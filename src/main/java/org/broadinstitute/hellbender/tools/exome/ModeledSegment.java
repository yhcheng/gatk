package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.Arrays;

/**
 * Class for storing a segment that was generated by a model (e.g., GATK CNV), though we do not necessarily know which model ahead of time.
 *
 * @author lichtens &lt;lichtens@broadinstitute.org&gt;
 */
public class ModeledSegment extends Segment<String> {

    public static final String NO_CALL = "";

    public ModeledSegment(final SimpleInterval interval, final String call, final long targetCount, final double segmentMeanInLogCR) {
        super(interval, targetCount, segmentMeanInLogCR, call);
        Utils.nonNull(interval, "The input interval cannot be null");
        Utils.nonNull(call, String.format("The input call cannot be null.  Use empty string, instead (\"%s\")", NO_CALL));
        ParamUtils.isFinite(segmentMeanInLogCR, "Segment Mean must be finite.");
    }

    public ModeledSegment(final SimpleInterval interval, final long targetCount, final double segmentMean) {
        this(interval, NO_CALL, targetCount, segmentMean);
    }

    /**
     *  Get segment mean in log2 space
     * @return
     */
    public double getSegmentMean() {
        return mean;
    }

    /**
     * Set the segment mean in log2 space
     *
     * @param segmentMean
     */
    public void setSegmentMean(final double segmentMean) {
        this.mean = segmentMean;
    }

    /**
     * Get the segment mean in non-logged space
     *
     * @return
     */
    public double getSegmentMeanInCRSpace() {
        return Math.pow(2, mean);
    }

    /**
     * Get the segment mean in logged space
     *
     * @return
     */
    public double getSegmentMeanInLog2CRSpace() {
        return getSegmentMean();
    }

    public void setSegmentMeanInCRSpace(final double segmentMeanInCRSpace) {
        this.mean = Math.log(segmentMeanInCRSpace)/Math.log(2);
    }

    @Override
    public String getContig() {return interval.getContig(); }

    @Override
    public int getStart() {return interval.getStart(); }

    @Override
    public int getEnd() {return interval.getEnd(); }

    public SimpleInterval getSimpleInterval() {
        return interval;
    }

    public void setSimpleInterval(final SimpleInterval simpleInterval) {
        this.interval = Utils.nonNull(simpleInterval, "The input interval cannot be null");
    }

    /**
     * Sets the call.
     */
    public void setCall(final String call) {
        this.call = Utils.nonNull(call, String.format("The input call cannot be null.  For no-calls use: \"%s\" ", NO_CALL));
    }

    public void setTargetCount(final long targetCount) {
        this.targetCount = ParamUtils.isPositiveOrZero(targetCount, "Number of targets must be positive or zero.");
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof ModeledSegment)) {
            return false;
        }

        final ModeledSegment modeledSegment = (ModeledSegment) o;
        return interval.equals(modeledSegment.interval) && call.equals(modeledSegment.call)
                && Math.abs(mean - modeledSegment.mean) < 2 * Math.ulp(mean)
                && targetCount == modeledSegment.targetCount;
    }

    @Override
    public int hashCode() {
        final int[] hashes = new int[]{ interval.hashCode(), call.hashCode(), Double.hashCode(mean),
                Long.hashCode(targetCount)};
        return Arrays.hashCode(hashes);
    }
}
