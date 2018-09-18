package org.broadinstitute.hellbender.utils.smithwaterman;


import com.github.yhcheng.ppc64native.smithwaterman.PPC64SmithWaterman;
//import com.ibm.ppc64native.smithwaterman.PPC64SmithWaterman;

import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWAlignerNativeBinding;
import org.broadinstitute.hellbender.exceptions.UserException;
/**
 * SmithWatermanIntelAligner class that converts instance of {@link SWAlignerNativeBinding} into a {@link SmithWatermanIntelAligner}
 * This is optimized for Intel Architectures and can fail if Machine does not support AVX and will throw {@link UserException}
 */

public final class SmithWatermanPPC64Aligner implements SmithWatermanAligner {



    private final SWAlignerNativeBinding smithWaterman = new PPC64SmithWaterman();

    /*
    * Generate SWAlignerWrapper instance
    */
    private final SWNativeAlignerWrapper alignerWrapper = new SWNativeAlignerWrapper(smithWaterman);


    /**
     * Create a new SW pairwise aligner, which is implementation of smith waterman aligner that's takes advantage of intel hardware optimizations.
     */

    public SmithWatermanPPC64Aligner() throws UserException.HardwareFeatureException {
        final boolean isSupported = smithWaterman.load(null);
        if (!isSupported) {
            throw new UserException.HardwareFeatureException("Machine does not support VSX SmithWaterman.");
        }
    }

    /**
     * Aligns the alternate sequence to the reference sequence
     *
     * @param reference  ref sequence
     * @param alternate  alt sequence
     */
    @Override
    public SmithWatermanAlignment align(final byte[] reference, final byte[] alternate, final SWParameters parameters, final SWOverhangStrategy overhangStrategy) {
//        logger.info("in SmithWatermanAlignment in ./gatk/src/main/java/org/broadinstitute/hellbender/utils/smithwaterman/SmithWatermanPPC64Aligner.java");
	return alignerWrapper.align(reference, alternate, parameters, overhangStrategy);
    }

    /**
     * Close the aligner
     */
    @Override
    public void close() {
        alignerWrapper.close();
    }
}
