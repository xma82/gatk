package org.broadinstitute.hellbender.utils.variant.writers;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.TreeRangeMap;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.List;

import static htsjdk.variant.vcf.VCFConstants.MAX_GENOTYPE_QUAL;

public class SomaticGVCFWriter {

    /** Where we'll ultimately write our VCF records */
    private final VariantContextWriter underlyingWriter;

    private final RangeMap<Double, Range<Double>> lodPartitions;

    /** fields updated on the fly during GVCFWriter operation */
    private int nextAvailableStart = -1;
    private String contigOfNextAvailableStart = null;
    private String sampleName = null;
    private TLODBlock currentBlock = null;

    /**
     * Create a new GVCF writer
     *
     * Should be a non-empty list of boundaries.  For example, suppose this variable is
     *
     * [A, B, C]
     *
     * We would partition our hom-ref sites into the following bands:
     *
     * X < A
     * A <= X < B
     * B <= X < C
     * X >= C
     *
     * @param underlyingWriter the ultimate destination of the GVCF records
     * @param lodPartitions     a list of GQ partitions, this list must be non-empty and every element must be larger than previous element
     */
    public SomaticGVCFWriter(final VariantContextWriter underlyingWriter, final List<Double> lodPartitions) {
        this.underlyingWriter = Utils.nonNull(underlyingWriter);
        this.lodPartitions = parsePartitions(lodPartitions);
    }

    /**
     * Create {@link HomRefBlock}s which will collectively accept variants of any genotype quality
     *
     * Each individual block covers a band of genotype qualities with the splits between bands occurring at values in {@code gqPartitions}.
     * There will be {@code gqPartitions.size() +1} bands produced covering the entire possible range of genotype qualities from 0 to {@link VCFConstants#MAX_GENOTYPE_QUAL}.
     *
     * @param lodPartitions proposed LOD partitions
     * @return a list of HomRefBlocks accepting bands of genotypes qualities split at the points specified in gqPartitions
     */
    @VisibleForTesting
    static RangeMap<Double,Range<Double>> parsePartitions(final List<Double> lodPartitions) {
        Utils.nonEmpty(lodPartitions);
        Utils.containsNoNull(lodPartitions, "The list of GQ partitions contains a null value");
        final RangeMap<Double, Range<Double>> result = TreeRangeMap.create();
        double lastThreshold = 0;
        for (final Double value : lodPartitions) {
            if (value < lastThreshold) {
                throw new IllegalArgumentException(String.format("The list of GQ partitions is out of order. Previous value is %d but the next is %d.", lastThreshold, value));
            } else if (value == lastThreshold) {
                throw new IllegalArgumentException(String.format("The value %d appears more than once in the list of GQ partitions.", value));
            }

            result.put(Range.closedOpen(lastThreshold, value), Range.closedOpen(lastThreshold, value));
            lastThreshold = value;
        }

        if (lastThreshold <= Double.POSITIVE_INFINITY) {
            result.put(Range.closedOpen(lastThreshold, Double.POSITIVE_INFINITY + 1), Range.closedOpen(lastThreshold,Double.POSITIVE_INFINITY + 1));
        }

        return result;
    }

    /**
     * Close this GVCF writer.  Finalizes any pending hom-ref blocks and emits those to the underlyingWriter as well
     */
    public void close() {
        try {
            emitCurrentBlock();
        } finally {
            underlyingWriter.close();
        }
    }

    public boolean checkError() {
        return underlyingWriter.checkError();
    }

    /**
     * Add hom-ref site from vc to this gVCF hom-ref state tracking, emitting any pending states if appropriate
     *
     * @param vc a non-null VariantContext
     * @param g  a non-null genotype from VariantContext
     * @return a VariantContext to be emitted, or null if non is appropriate
     */
    protected VariantContext addHomRefSite(final VariantContext vc, final Genotype g) {

        if (nextAvailableStart != -1) {
            // don't create blocks while the hom-ref site falls before nextAvailableStart (for deletions)
            if (vc.getStart() <= nextAvailableStart && vc.getContig().equals(contigOfNextAvailableStart)) {
                return null;
            }
            // otherwise, reset to non-relevant
            nextAvailableStart = -1;
            contigOfNextAvailableStart = null;
        }

        final VariantContext result;
        if (genotypeCanBeMergedInCurrentBlock(g)) {
            currentBlock.add(vc.getStart(), g);
            result = null;
        } else {
            result = currentBlock != null ? currentBlock.toVariantContext(sampleName): null;
            currentBlock = createNewBlock(vc, g);
        }
        return result;
    }

    private boolean genotypeCanBeMergedInCurrentBlock(final Genotype g) {
        return currentBlock != null
                && currentBlock.withinBounds(g.getAttributeAsDouble(GATKVCFConstants.TUMOR_LOD_KEY, Double.NEGATIVE_INFINITY));
    }

    /**
     * Flush the current hom-ref block, if necessary, to the underlying writer, and reset the currentBlock to null
     */
    private void emitCurrentBlock() {
        if (currentBlock != null) {
            underlyingWriter.add(currentBlock.toVariantContext(sampleName));
            currentBlock = null;
        }
    }


    /**
     * Helper function to create a new HomRefBlock from a variant context and current genotype
     *
     * @param vc the VariantContext at the site where want to start the band
     * @param g  the genotype of the sample from vc that should be used to initialize the block
     * @return a newly allocated and initialized block containing g already
     */
    private TLODBlock createNewBlock(final VariantContext vc, final Genotype g) {
        // figure out the GQ limits to use based on the GQ of g
        final double lod = g.getAttributeAsDouble(GATKVCFConstants.TUMOR_LOD_KEY, Double.NEGATIVE_INFINITY);
        final Range<Double> partition = lodPartitions.get(lod);

        if( partition == null) {
            throw new GATKException("GQ " + g + " from " + vc + " didn't fit into any partition");
        }

        // create the block, add g to it, and return it for use
        final TLODBlock block = new TLODBlock(vc, partition.lowerEndpoint(), partition.upperEndpoint());
        block.add(vc.getStart(), g);
        return block;
    }

    /**
     * Add a VariantContext to this writer for emission
     *
     * Requires that the VC have exactly one genotype
     *
     * @param vc a non-null VariantContext
     */
    public void add(VariantContext vc) {
        Utils.nonNull(vc);
        Utils.validateArg(vc.hasGenotypes(), "GVCF assumes that the VariantContext has genotypes");
        Utils.validateArg(vc.getGenotypes().size() == 1, () -> "GVCF assumes that the VariantContext has exactly one genotype but saw " + vc.getGenotypes().size());

        if (sampleName == null) {
            sampleName = vc.getGenotype(0).getSampleName();
        }

        if (currentBlock != null && !currentBlock.isContiguous(vc)) {
            // we've made a non-contiguous step (across interval, onto another chr), so finalize
            emitCurrentBlock();
        }

        final Genotype g = vc.getGenotype(0);
        if (g.isHomRef() && vc.hasAlternateAllele(Allele.NON_REF_ALLELE) && vc.isBiallelic()) {
            // create bands
            final VariantContext maybeCompletedBand = addHomRefSite(vc, g);
            if (maybeCompletedBand != null) {
                underlyingWriter.add(maybeCompletedBand);
            }
        } else {
            // g is variant, so flush the bands and emit vc
            emitCurrentBlock();
            nextAvailableStart = vc.getEnd();
            contigOfNextAvailableStart = vc.getContig();
            underlyingWriter.add(vc);
        }

    }
}
