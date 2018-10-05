package org.broadinstitute.hellbender.utils.variant.writers;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;

/**
 * Helper class for calculating a somatic LOD band in the SomaticGVCF writer
 *
 * A band contains LOD and DP values for a contiguous stretch of hom-ref genotypes,
 * and provides summary information about the entire block of genotypes.
 *
 * Genotypes within the TLODBlock are restricted to hom-ref genotypes within a band of LOD scores
 */
final class TLODBlock implements Locatable {

    private final VariantContext startingVC;
    private final double minLOD, maxLOD;
    private final List<Integer> DPs = new ArrayList<>();
    private final Allele ref;

    private int end;
    private double maxBlockLOD = Double.NEGATIVE_INFINITY;
    //given that we're using the same LOD calculation as for
    // variants, more confident reference sites will have lower LOD so the most conservative thing is to take the max LOD of the block

    /**
     * Create a new HomRefBlock
     *
     * @param startingVC the VariantContext that starts this band (for starting position information)
     * @param lowerLODBound the lowerLODBound (inclusive) to use in this band
     * @param upperLODBound the upperLODBound (exclusive) to use in this band
     */
    public TLODBlock(final VariantContext startingVC, final double lowerLODBound, final double upperLODBound) {
        Utils.nonNull(startingVC, "startingVC cannot be null");
        if ( lowerLODBound > upperLODBound ) { throw new IllegalArgumentException("bad lowerLODBound " + lowerLODBound + " as it's >= upperLODBound " + upperLODBound); }

        this.startingVC = startingVC;
        this.end = getStart() - 1;
        this.ref = startingVC.getReference();
        this.minLOD = lowerLODBound;
        this.maxLOD = upperLODBound;
    }

    /**
     * Convert a HomRefBlock into a VariantContext
     *
     * @param sampleName sample name to give this variant context
     * @return a VariantContext representing the gVCF encoding for this block.
     * It will return {@code null} if input {@code block} is {@code null}, indicating that there
     * is no variant-context to be output into the VCF.
     */
    public VariantContext toVariantContext(String sampleName) {
        final VariantContextBuilder vcb = new VariantContextBuilder(getStartingVC());
        vcb.attributes(new LinkedHashMap<>(2)); // clear the attributes
        vcb.stop(getEnd());
        vcb.attribute(VCFConstants.END_KEY, getEnd());
        final Genotype genotype = createHomRefGenotype(sampleName);

        return vcb.genotypes(genotype).make();
    }

    // create a single Genotype with GQ and DP annotations
    private Genotype createHomRefGenotype(String sampleName) {
        final GenotypeBuilder gb = new GenotypeBuilder(sampleName, Collections.nCopies(2, getRef()));  //FIXME: for somatic stuff we output the genotype as diploid because that's familiar for human
        gb.noAD().noPL().noAttributes(); // clear all attributes

        gb.attribute(GATKVCFConstants.TUMOR_LOD_KEY, maxBlockLOD);
        gb.DP(getMedianDP());
        gb.attribute(GATKVCFConstants.MIN_DP_FORMAT_KEY, getMinDP());

        return gb.make();
    }

    /**
     * Add information from this Genotype to this band.
     *
     * @param pos Current genomic position. Must be 1 base after the previous position
     * @param genotype A non-null Genotype with GQ and DP attributes
     */
    public void add(final int pos, final Genotype genotype) {
        Utils.nonNull(genotype, "genotype cannot be null");
        if ( pos != end + 1 ) { throw new IllegalArgumentException("adding genotype at pos " + pos + " isn't contiguous with previous end " + end); }
        // Make sure the LOD is within the bounds of this band
        if ( !withinBounds(genotype.getAttributeAsDouble(GATKVCFConstants.TUMOR_LOD_KEY, Double.NEGATIVE_INFINITY))) {
            throw new IllegalArgumentException("cannot add a genotype with LOD=" + genotype.getGQ() + " because it's not within bounds ["
                    + this.getLODLowerBound() + ',' + this.getLODUpperBound() + ')');
        }

        double currentLOD = genotype.getAttributeAsDouble(GATKVCFConstants.TUMOR_LOD_KEY, Double.NEGATIVE_INFINITY);
        if( maxBlockLOD == Double.NEGATIVE_INFINITY || currentLOD > maxBlockLOD) {
            maxBlockLOD = currentLOD;
        }

        end = pos;
        DPs.add(Math.max(genotype.getDP(), 0)); // DP must be >= 0
    }

    /**
     * Is the LOD value within the bounds of this LOD (LOD >= minLOD && LOD < maxLOD)
     * @param LOD the LOD value to test
     * @return true if within bounds, false otherwise
     */
    public boolean withinBounds(final double LOD) {
        return LOD >= minLOD && LOD < maxLOD;
    }

    /** Get the min DP observed within this band */
    public int getMinDP() {
        return Collections.min(DPs);
    }

    /** Get the median DP observed within this band
     * If there are an even number of DPs recorded in this band the median is the mean of the two middle values */
    public int getMedianDP() {
        return (int) Math.round(MathUtils.median(DPs));
    }

    /** Get the min PLs observed within this band, can be null if no PLs have yet been observed */
    public double getMaxBlockLOD() {
        return maxBlockLOD;
    }

    double getLODUpperBound() {
        return maxLOD;
    }
    double getLODLowerBound() {
        return minLOD;
    }

    public boolean isContiguous(final VariantContext vc) {
        return (vc.getEnd() == getEnd() + 1) && startingVC.getContig().equals(vc.getContig());
    }

    public VariantContext getStartingVC() {
        return startingVC;
    }

    @Override
    public String getContig() {
        return startingVC.getContig();
    }

    @Override
    public int getStart() {
        return startingVC.getStart();
    }

    @Override
    public int getEnd() {
        return end;
    }

    public Allele getRef() {
        return ref;
    }

    public int getSize() {
        return getEnd() - getStart() + 1;
    }

    @Override
    public String toString() {
        return "HomRefBlock{" +
                "minLOD=" + minLOD +
                ", maxLOD=" + maxLOD +
                '}';
    }
}

