package org.broadinstitute.hellbender.tools.walkers.mutect;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.tools.walkers.genotyper.PloidyModel;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceModel;
import org.broadinstitute.hellbender.tools.walkers.variantutils.PosteriorProbabilitiesUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;
import java.util.stream.DoubleStream;

import static org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceModel.isAltAfterAssembly;
import static org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceModel.isAltBeforeAssembly;
import static org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils.getOverlappingVariantContext;

public class SomaticReferenceConfidenceModel {

    private final SampleList samples;
    private final int indelInformativeDepthIndelSize;
    private final SomaticGenotypingEngine genotypingEngine;

    /**
     * Holds information about a genotype call of a single sample reference vs. any non-ref event
     *
     * IMPORTANT PERFORMANCE NOTE!!! Allowing direct field access (within this class only) speeds up
     * the HaplotypeCaller by ~10% vs. accessing the fields indirectly via setters, as seen in a profiler.
     */
    @VisibleForTesting
    public static final class SomaticRefVsAnyResult {
        /**
         * The genotype likelihoods for ref/ref ref/non-ref non-ref/non-ref
         *
         * Fields are visible because direct field access for this particular class has a major performance
         * impact on the HaplotypeCaller, as noted above, and the class itself is nested within
         * ReferenceConfidenceModel anyway.
         */
        PerAlleleCollection<Double> lods;

        int refDepth = 0;
        int nonRefDepth = 0;

        /**
         * Creates a new ref-vs-alt result indicating the genotype likelihood vector capacity.
         */
        public SomaticRefVsAnyResult() {
            lods = new PerAlleleCollection<>(PerAlleleCollection.Type.REF_AND_ALT);
        }

        /**
         * @return Get the DP (sum of AD values)
         */
        int getDP() {
            return refDepth + nonRefDepth;
        }

        /**
         * Return the AD fields. Returns a newly allocated array every time.
         */
        int[] getAD() {
            return new int[]{refDepth, nonRefDepth};
        }
    }

    /**
     * Create a new ReferenceConfidenceModel
     *
     * @param samples the list of all samples we'll be considering with this model
     * @param header the SAMFileHeader describing the read information (used for debugging)
     * @param indelInformativeDepthIndelSize the max size of indels to consider when calculating indel informative depths
     */
    public SomaticReferenceConfidenceModel(final SampleList samples,
                                    final SAMFileHeader header,
                                    final int indelInformativeDepthIndelSize,
                                    final SomaticGenotypingEngine genotypingEngine){
        this.samples = samples;
        this.indelInformativeDepthIndelSize = indelInformativeDepthIndelSize;
        this.genotypingEngine = genotypingEngine;
    }

    /**
     * Calculate the genotype likelihoods for the sample in pileup for being hom-ref contrasted with being ref vs. alt
     *
     * @param ploidy target sample ploidy.
     * @param pileup the read backed pileup containing the data we want to evaluate
     * @param refBase the reference base at this pileup position
     * @param qual the min base quality for a read in the pileup at the pileup position to be included in the calculation
     * @param hqSoftClips running average data structure (can be null) to collect information about the number of high quality soft clips
     * @return a RefVsAnyResult genotype call.
     */
    public SomaticRefVsAnyResult calcGenotypeLikelihoodsOfRefVsAny(final int ploidy,
                                                            final ReadPileup pileup,
                                                            final byte refBase,
                                                            final byte qual,
                                                            final MathUtils.RunningAverage hqSoftClips,
                                                            final boolean readsWereRealigned) {

        final SomaticRefVsAnyResult result = new SomaticRefVsAnyResult();
        Map<String, List<GATKRead>> perSampleReadMap = new HashMap<>();
        perSampleReadMap.put(samples.getSample(0), pileup.getReads());
        ReadLikelihoods readLikelihoods = new ReadLikelihoods(samples, new IndexedAlleleList(Arrays.asList(Allele.create(refBase,true), Allele.NON_REF_ALLELE)), perSampleReadMap);
        final Iterator<PileupElement> pileupIter = pileup.iterator();
        for (int i = 0; i < pileup.size(); i++) {
            final PileupElement element = pileupIter.next();
            final boolean isAlt = readsWereRealigned ? isAltAfterAssembly(element, refBase) : isAltBeforeAssembly(element, refBase);
            final double referenceLikelihood;
            final double nonRefLikelihood;
            if (isAlt) {
                nonRefLikelihood = QualityUtils.qualToProbLog10(qual);
                referenceLikelihood = QualityUtils.qualToErrorProbLog10(qual) + MathUtils.LOG10_ONE_THIRD;
                result.nonRefDepth++;
            } else {
                referenceLikelihood = QualityUtils.qualToProbLog10(qual);
                nonRefLikelihood = QualityUtils.qualToErrorProbLog10(qual) + MathUtils.LOG10_ONE_THIRD;
                result.refDepth++;
            }
            readLikelihoods.sampleMatrix(0).set(0, i, nonRefLikelihood);
        }
        result.lods = genotypingEngine.somaticLog10Odds(readLikelihoods.sampleMatrix(0));
        return result;
    }

    public List<VariantContext> calculateRefConfidence(final Haplotype refHaplotype,
                                                       final Collection<Haplotype> calledHaplotypes,
                                                       final SimpleInterval paddedReferenceLoc,
                                                       final AssemblyRegion activeRegion,
                                                       final ReadLikelihoods<Haplotype> readLikelihoods,
                                                       final PloidyModel ploidyModel,
                                                       final List<VariantContext> variantCalls) {
        return calculateRefConfidence(refHaplotype, calledHaplotypes, paddedReferenceLoc, activeRegion, readLikelihoods,
                ploidyModel, variantCalls, false, Collections.emptyList());
    }

    /**
     * Calculate the reference confidence for a single sample given the its read data
     *
     * Returns a list of variant contexts, one for each position in the {@code activeRegion.getLoc()}, each containing
     * detailed information about the certainty that the sample is hom-ref for each base in the region.
     *
     *
     *
     * @param refHaplotype the reference haplotype, used to get the reference bases across activeRegion.getLoc()
     * @param calledHaplotypes a list of haplotypes that segregate in this region, for realignment of the reads in the
     *                         readLikelihoods, corresponding to each reads best haplotype.  Must contain the refHaplotype.
     * @param paddedReferenceLoc the location of refHaplotype (which might be larger than activeRegion.getLoc())
     * @param activeRegion the active region we want to get the reference confidence over
     * @param readLikelihoods a map from a single sample to its PerReadAlleleLikelihoodMap for each haplotype in calledHaplotypes
     * @param ploidyModel indicate the ploidy of each sample in {@code stratifiedReadMap}.
     * @param variantCalls calls made in this region.  The return result will contain any variant call in this list in the
     *                     correct order by genomic position, and any variant in this list will stop us emitting a ref confidence
     *                     under any position it covers (for snps and insertions that is 1 bp, but for deletions its the entire ref span)
     * @return an ordered list of variant contexts that spans activeRegion.getLoc() and includes both reference confidence
     *         contexts as well as calls from variantCalls if any were provided
     */
    public List<VariantContext> calculateRefConfidence(final Haplotype refHaplotype,
                                                       final Collection<Haplotype> calledHaplotypes,
                                                       final SimpleInterval paddedReferenceLoc,
                                                       final AssemblyRegion activeRegion,
                                                       final ReadLikelihoods<Haplotype> readLikelihoods,
                                                       final PloidyModel ploidyModel,
                                                       final List<VariantContext> variantCalls,
                                                       final boolean applyPriors,
                                                       final List<VariantContext> VCpriors) {
        Utils.nonNull(refHaplotype, "refHaplotype cannot be null");
        Utils.nonNull(calledHaplotypes, "calledHaplotypes cannot be null");
        Utils.validateArg(calledHaplotypes.contains(refHaplotype), "calledHaplotypes must contain the refHaplotype");
        Utils.nonNull(paddedReferenceLoc, "paddedReferenceLoc cannot be null");
        Utils.nonNull(activeRegion, "activeRegion cannot be null");
        Utils.nonNull(readLikelihoods, "readLikelihoods cannot be null");
        Utils.validateArg(readLikelihoods.numberOfSamples() == 1, () -> "readLikelihoods must contain exactly one sample but it contained " + readLikelihoods.numberOfSamples());
        Utils.validateArg( refHaplotype.length() == activeRegion.getExtendedSpan().size(), () -> "refHaplotype " + refHaplotype.length() + " and activeRegion location size " + activeRegion.getSpan().size() + " are different");
        Utils.nonNull(ploidyModel, "the ploidy model cannot be null");
        final int ploidy = ploidyModel.samplePloidy(0); // the first sample = the only sample in reference-confidence mode.

        final SimpleInterval refSpan = activeRegion.getSpan();
        final List<ReadPileup> refPileups = AssemblyBasedCallerUtils.getPileupsOverReference(refHaplotype, calledHaplotypes, paddedReferenceLoc, activeRegion, refSpan, readLikelihoods, samples);
        final byte[] ref = refHaplotype.getBases();
        final List<VariantContext> results = new ArrayList<>(refSpan.size());
        final String sampleName = readLikelihoods.getSample(0);

        final int globalRefOffset = refSpan.getStart() - activeRegion.getExtendedSpan().getStart();
        for ( final ReadPileup pileup : refPileups ) {
            final Locatable curPos = pileup.getLocation();
            final int offset = curPos.getStart() - refSpan.getStart();

            final VariantContext overlappingSite = getOverlappingVariantContext(curPos, variantCalls);
            if ( overlappingSite != null && overlappingSite.getStart() == curPos.getStart() ) {
                    results.add(overlappingSite);
            } else {
                // otherwise emit a reference confidence variant context
                results.add(makeReferenceConfidenceVariantContext(ploidy, ref, sampleName, globalRefOffset, pileup, curPos, offset, readLikelihoods));
            }
        }

        return results;
    }

    public VariantContext makeReferenceConfidenceVariantContext(final int ploidy,
                                                                final byte[] ref,
                                                                final String sampleName,
                                                                final int globalRefOffset,
                                                                final ReadPileup pileup,
                                                                final Locatable curPos,
                                                                final int offset,
                                                                final ReadLikelihoods readLikelihoods) {
        // Assume infinite population on a single sample.
        final int refOffset = offset + globalRefOffset;
        final byte refBase = ref[refOffset];
        final SomaticRefVsAnyResult result = calcGenotypeLikelihoodsOfRefVsAny(ploidy, pileup, refBase, (byte)6, null, true);

        final Allele refAllele = Allele.create(refBase, true);
        final List<Allele> refSiteAlleles = Arrays.asList(refAllele, Allele.NON_REF_ALLELE);
        final VariantContextBuilder vcb = new VariantContextBuilder("HC", curPos.getContig(), curPos.getStart(), curPos.getStart(), refSiteAlleles);
        final GenotypeBuilder gb = new GenotypeBuilder(sampleName, GATKVariantContextUtils.homozygousAlleleList(refAllele, ploidy));
        gb.AD(result.getAD());
        gb.DP(result.getDP());

        // genotype likelihood calculation
        //final int nIndelInformativeReads = calcNIndelInformativeReads(pileup, refOffset, ref, indelInformativeDepthIndelSize);
        //final Double indelLod =

        // now that we have the SNP and indel GLs, we take the one with the least confidence,
        // as this is the most conservative estimate of our certainty that we are hom-ref.
        // For example, if the SNP PLs are 0,10,100 and the indel PLs are 0,100,1000
        // we are very certain that there's no indel here, but the SNP confidence imply that we are
        // far less confident that the ref base is actually the only thing here.  So we take 0,10,100
        // as our GLs for the site.
        gb.attribute(GATKVCFConstants.TUMOR_LOD_KEY, result.lods.get(Allele.NON_REF_ALLELE));


        return vcb.genotypes(gb.make()).make();
    }

}
