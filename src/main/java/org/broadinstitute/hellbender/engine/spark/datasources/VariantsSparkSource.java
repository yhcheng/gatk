package org.broadinstitute.hellbender.engine.spark.datasources;

import htsjdk.samtools.seekablestream.SeekableStreamFactory;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.io.LongWritable;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;
import org.broadinstitute.hellbender.utils.variant.VariantContextVariantAdapter;
import org.seqdoop.hadoop_bam.VCFInputFormat;
import org.seqdoop.hadoop_bam.VariantContextWritable;
import org.seqdoop.hadoop_bam.util.BGZFCodec;
import org.seqdoop.hadoop_bam.util.BGZFEnhancedGzipCodec;
import org.seqdoop.hadoop_bam.util.VCFHeaderReader;

import java.io.IOException;
import java.util.List;

/**
 * VariantsSparkSource loads Variants from files serially (using FeatureDataSource<VariantContext>) or in parallel
 * using Hadoop-BAM.
 */
public final class VariantsSparkSource {
    private final JavaSparkContext ctx;

    public VariantsSparkSource(JavaSparkContext ctx) {
        this.ctx = ctx;
    }

    /**
     * Loads variants in parallel using Hadoop-BAM for vcfs and bcfs.
     * @param vcf file to load variants from.
     * @param intervals intervals of variants to include, or null if all should be included.
     * @return JavaRDD<GATKVariant> of variants from the variants file specified in vcf.
     */
    public JavaRDD<GATKVariant> getParallelVariants(final String vcf, final List<SimpleInterval> intervals) {
        return getParallelVariantContexts(vcf, intervals)
                .filter(vc -> vc.getCommonInfo() != null)
                .map(vc -> VariantContextVariantAdapter.sparkVariantAdapter(vc));
    }

    /**
     * Loads variants from multiple inputs in parallel using Hadoop-BAM for vcfs and bcfs.
     * @param vcfs List of input files to load variants from.
     * @param intervals intervals of variants to include, or null if all should be included.
     * @return JavaRDD<VariantContext> rdd containing the union of all variants from the variant
     * files specified in vcfs.
     */
    public JavaRDD<GATKVariant> getParallelVariants(final List<String> vcfs, final List<SimpleInterval> intervals) {
        return vcfs.parallelStream()
                .map(vcf -> getParallelVariants(vcf, intervals))
                .reduce(ctx.emptyRDD(), (result, rdd) -> result.union(rdd)
        );
    }

    /**
     * Loads variants in parallel using Hadoop-BAM for vcfs and bcfs.
     * @param vcf file to load variants from.
     * @param intervals intervals of variants to include, or null if all should be included.
     * @return JavaRDD<VariantContext> of variants from all files.
     */
    public JavaRDD<VariantContext> getParallelVariantContexts(final String vcf, final List<SimpleInterval> intervals) {
        Configuration conf = new Configuration();
        conf.setStrings("io.compression.codecs", BGZFEnhancedGzipCodec.class.getCanonicalName(),
                BGZFCodec.class.getCanonicalName());
        if (intervals != null && !intervals.isEmpty()) {
            VCFInputFormat.setIntervals(conf, intervals);
        }
        final JavaPairRDD<LongWritable, VariantContextWritable> rdd2 = ctx.newAPIHadoopFile(
                vcf, VCFInputFormat.class, LongWritable.class, VariantContextWritable.class,
                conf);
        return rdd2.map(v1 -> v1._2().get());
    }

    public static VCFHeader getHeader(String filePath) {
        try {
            return VCFHeaderReader.readHeaderFrom(SeekableStreamFactory.getInstance().getStreamFor(filePath));
        } catch (IOException e) {
            throw new UserException("Failed to read VCF header from " + filePath + "\n Caused by:" + e.getMessage(), e);
        }
    }
}
