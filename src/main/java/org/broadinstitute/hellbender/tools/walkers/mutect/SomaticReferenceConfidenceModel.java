package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceModel;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

public class SomaticReferenceConfidenceModel extends ReferenceConfidenceModel {
    @Override
    /**
     * Calculate the genotype likelihoods for the sample in pileup for being hom-ref contrasted with being ref vs. alt
     *
     * @param ploidy target sample ploidy.
     * @param pileup the read backed pileup containing the data we want to evaluate
     * @param refBase the reference base at this pileup position
     * @param minBaseQual the min base quality for a read in the pileup at the pileup position to be included in the calculation
     * @param hqSoftClips running average data structure (can be null) to collect information about the number of high quality soft clips
     * @return a RefVsAnyResult genotype call.
     */
    public RefVsAnyResult calcGenotypeLikelihoodsOfRefVsAny(final int ploidy,
                                                            final ReadPileup pileup,
                                                            final byte refBase,
                                                            final byte minBaseQual,
                                                            final MathUtils.RunningAverage hqSoftClips,
                                                            final boolean readsWereRealigned) {
}
