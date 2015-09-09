package org.broadinstitute.hellbender.utils;

import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import htsjdk.variant.vcf.VCFHeader;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;


/**
 * SampleUtils is a static class (no instantiation allowed!) with some utility methods for getting samples
 * quality scores.
 *
 * @author ebanks
 */
public class SampleUtils {
    /**
     * Private constructor.  No instantiating this class!
     */
    private SampleUtils() {}

    public static Set<String> getSampleList(Map<String, VCFHeader> headers, GATKVariantContextUtils.GenotypeMergeType mergeOption) {
        Set<String> samples = new TreeSet<String>();
        for ( Map.Entry<String, VCFHeader> val : headers.entrySet() ) {
            VCFHeader header = val.getValue();
            for ( String sample : header.getGenotypeSamples() ) {
                samples.add(GATKVariantContextUtils.mergedSampleName(val.getKey(), sample, mergeOption == GATKVariantContextUtils.GenotypeMergeType.UNIQUIFY));
            }
        }

        return samples;
    }
}
