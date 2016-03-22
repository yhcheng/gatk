package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfile;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.activityprofile.BandPassActivityProfile;
import org.broadinstitute.hellbender.utils.downsampling.DownsamplingMethod;
import org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;


/**
 * Call germline SNPs and indels via local re-assembly of haplotypes
 *
 * <p>The HaplotypeCaller is capable of calling SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region. In other words, whenever the program encounters a region showing signs of variation, it discards the existing mapping information and completely reassembles the reads in that region. This allows the HaplotypeCaller to be more accurate when calling regions that are traditionally difficult to call, for example when they contain different types of variants close to each other. It also makes the HaplotypeCaller much better at calling indels than position-based callers like UnifiedGenotyper.</p>
 *
 * <p>In the so-called GVCF mode used for scalable variant calling in DNA sequence data, HaplotypeCaller runs per-sample to generate an intermediate genomic gVCF (gVCF), which can then be used for joint genotyping of multiple samples in a very efficient way, which enables rapid incremental processing of samples as they roll off the sequencer, as well as scaling to very large cohort sizes (e.g. the 92K exomes of ExAC).</p>
 *
 * <p>In addition, HaplotypeCaller is able to handle non-diploid organisms as well as pooled experiment data. Note however that the algorithms used to calculate variant likelihoods is not well suited to extreme allele frequencies (relative to ploidy) so its use is not recommended for somatic (cancer) variant discovery. For that purpose, use MuTect2 instead.</p>
 *
 * <p>Finally, HaplotypeCaller is also able to correctly handle the splice junctions that make RNAseq a challenge for most variant callers.</p>
 *
 * <h3>How HaplotypeCaller works</h3>
 *
 * <br />
 * <h4>1. Define active regions </h4>
 *
 * <p>The program determines which regions of the genome it needs to operate on, based on the presence of significant
 * evidence for variation.</p>
 *
 * <br />
 * <h4>2. Determine haplotypes by assembly of the active region </h4>
 *
 * <p>For each ActiveRegion, the program builds a De Bruijn-like graph to reassemble the ActiveRegion, and identifies
 * what are the possible haplotypes present in the data. The program then realigns each haplotype against the reference
 * haplotype using the Smith-Waterman algorithm in order to identify potentially variant sites. </p>
 *
 * <br />
 * <h4>3. Determine likelihoods of the haplotypes given the read data </h4>
 *
 * <p>For each ActiveRegion, the program performs a pairwise alignment of each read against each haplotype using the
 * PairHMM algorithm. This produces a matrix of likelihoods of haplotypes given the read data. These likelihoods are
 * then marginalized to obtain the likelihoods of alleles for each potentially variant site given the read data.   </p>
 *
 * <br />
 * <h4>4. Assign sample genotypes </h4>
 *
 * <p>For each potentially variant site, the program applies Bayes' rule, using the likelihoods of alleles given the
 * read data to calculate the likelihoods of each genotype per sample given the read data observed for that
 * sample. The most likely genotype is then assigned to the sample.    </p>
 *
 * <h3>Input</h3>
 * <p>
 * Input bam file(s) from which to make calls
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * Either a VCF or gVCF file with raw, unfiltered SNP and indel calls. Regular VCFs must be filtered either by variant
 * recalibration (best) or hard-filtering before use in downstream analyses. If using the reference-confidence model
 * workflow for cohort analysis, the output is a GVCF file that must first be run through GenotypeGVCFs and then
 * filtering before further analysis.
 * </p>
 *
 * <h3>Usage examples</h3>
 *
 * <p>These are example commands that show how to run HaplotypeCaller for typical use cases. Square brackets ("[ ]")
 * indicate optional arguments. Note that parameter values shown here may not be the latest recommended; see the
 * Best Practices documentation for detailed recommendations. </p>
 *
 * <br />
 * <h4>Single-sample GVCF calling on DNAseq (for `-ERC GVCF` cohort analysis workflow)</h4>
 * <pre>
 *   java -jar GenomeAnalysisTK.jar \
 *     -R reference.fasta \
 *     -T HaplotypeCaller \
 *     -I sample1.bam \
 *     --emitRefConfidence GVCF \
 *     [--dbsnp dbSNP.vcf] \
 *     [-L targets.interval_list] \
 *     -o output.raw.snps.indels.g.vcf
 * </pre>
 *
 * <h4>Single-sample GVCF calling on DNAseq with allele-specific annotations (for allele-specific cohort analysis workflow)</h4>
 * <pre>
 *   java -jar GenomeAnalysisTK.jar \
 *     -R reference.fasta \
 *     -T HaplotypeCaller \
 *     -I sample1.bam \
 *     --emitRefConfidence GVCF \
 *     [--dbsnp dbSNP.vcf] \
 *     [-L targets.interval_list] \
 *     -G Standard -G AS_Standard \
 *     -o output.raw.snps.indels.AS.g.vcf
 * </pre>
 *
 * <h4>Variant-only calling on DNAseq</h4>
 * <pre>
 *   java -jar GenomeAnalysisTK.jar \
 *     -R reference.fasta \
 *     -T HaplotypeCaller \
 *     -I sample1.bam [-I sample2.bam ...] \
 *     [--dbsnp dbSNP.vcf] \
 *     [-stand_call_conf 30] \
 *     [-stand_emit_conf 10] \
 *     [-L targets.interval_list] \
 *     -o output.raw.snps.indels.vcf
 * </pre>
 *
 * <h4>Variant-only calling on RNAseq</h4>
 * <pre>
 *   java -jar GenomeAnalysisTK.jar \
 *     -R reference.fasta \
 *     -T HaplotypeCaller \
 *     -I sample1.bam \
 *     [--dbsnp dbSNP.vcf] \
 *     -stand_call_conf 20 \
 *     -stand_emit_conf 20 \
 *     -o output.raw.snps.indels.vcf
 * </pre>
 *
 * <h3>Caveats</h3>
 * <ul>
 * <li>We have not yet fully tested the interaction between the GVCF-based calling or the multisample calling and the
 * RNAseq-specific functionalities. Use those in combination at your own risk.</li>
 * <li>Many users have reported issues running HaplotypeCaller with the -nct argument, so we recommend using Queue to
 * parallelize HaplotypeCaller instead of multithreading.</li>
 * </ul>
 *
 * <h3>Special note on ploidy</h3>
 * <p>This tool is able to handle almost any ploidy (except very high ploidies in large pooled experiments); the ploidy can be specified using the -ploidy argument for non-diploid organisms.</p>
 *
 * <h3>Additional Notes</h3>
 * <ul>
 *     <li>When working with PCR-free data, be sure to set `-pcr_indel_model NONE` (see argument below).</li>
 *     <li>When running in `-ERC GVCF` or `-ERC BP_RESOLUTION` modes, the emitting and calling confidence thresholds
 *     are automatically set to 0. This cannot be overridden by the command line. The thresholds can be set manually
 *     to the desired levels in the next step of the workflow (GenotypeGVCFs)</li>
 * </ul>
 *
 */
@CommandLineProgramProperties(
        summary = "Call germline SNPs and indels via local re-assembly of haplotypes",
        oneLineSummary = "Call germline SNPs and indels via local re-assembly of haplotypes",
        programGroup = VariantProgramGroup.class
)
public class HaplotypeCaller extends ReadWindowWalker {

    @ArgumentCollection
    private HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();

    /**
     * A raw, unfiltered, highly sensitive callset in VCF format.
     */
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "File to which variants should be written")
    public String outputVCF = null;

    @Argument(fullName = "activeRegionsOnly", shortName = "AR", doc = "Call only active regions", optional = true)
    public boolean activeRegionsOnly = false;

    private VariantContextWriter vcfWriter;

    private HaplotypeCallerEngine hcEngine;

    @Override
    public boolean requiresReference() { return true; }

    @Override
    protected int defaultWindowSize() { return activeRegionsOnly ? 5000 : 300; }

    @Override
    protected int defaultWindowPadding() { return 100; }

    @Override
    public CountingReadFilter makeReadFilter() {
        return HaplotypeCallerEngine.makeStandardHCReadFilter(hcArgs, getHeaderForReads());
    }

    @Override
    public void onTraversalStart() {
        hcEngine = new HaplotypeCallerEngine(hcArgs, getHeaderForReads(), referenceArguments.getReferenceFileName());

        // The HC engine will make the right kind (VCF or GVCF) of writer for us
        vcfWriter = hcEngine.makeVCFWriter(outputVCF, getHeaderForReads().getSequenceDictionary());
        hcEngine.writeHeader(vcfWriter);
    }

    @Override
    public void apply( ReadWindow window, ReferenceContext referenceContext, FeatureContext featureContext ) {
        Iterable<AssemblyRegion> assemblyRegions = activeRegionsOnly ? readWindowToActiveRegions(window, referenceContext) :
                                                                       readWindowToAssemblyRegions(window);

        for ( final AssemblyRegion assemblyRegion : assemblyRegions ) {
            logger.debug("Processing assembly region at " + assemblyRegion.getSpan() + " isActive: " + assemblyRegion.isActive() + " numReads: " + assemblyRegion.getReads().size());

            List<VariantContext> callsInRegion = hcEngine.callRegion(assemblyRegion, featureContext);

            for ( final VariantContext call : callsInRegion ) {
                // Only include calls that start within the window (as opposed to the padded regions around it)
                if ( window.startsWithinWindow(call) ) {
                    vcfWriter.add(call);
                }
            }
        }
    }

    /**
     * Convert a {@link ReadWindow} into a List of {@link AssemblyRegion} for processing by the HaplotypeCaller,
     * covering the span of the original window
     *
     * @param window window to convert into AssemblyRegions
     * @return a List of {@link AssemblyRegion} for processing by the HaplotypeCaller, covering the span of
     *         the original window
     */
    private Iterable<AssemblyRegion> readWindowToAssemblyRegions( final ReadWindow window ) {
        List<AssemblyRegion> regions = new ArrayList<>();

        // Currently, we convert each window 1:1 into an equivalent AssemblyRegion, marked as active.
        // We could call into {@link HaplotypeCallerEngine#isActive} here (or similar) to divide each
        // window into active and inactive regions, but we have to prove that the cost of doing so is
        // worth the gains.
        final AssemblyRegion monolithicRegion = new AssemblyRegion(window.getInterval(), windowPadding, getHeaderForReads());
        for ( final GATKRead read : window ) {
            monolithicRegion.add(read);
        }

        regions.add(monolithicRegion);
        return regions;
    }

    private Iterable<AssemblyRegion> readWindowToActiveRegions( final ReadWindow window, final ReferenceContext referenceContext  ) {
        final List<GATKRead> windowReads = new ArrayList<>();
        for ( final GATKRead read : window ) {
            windowReads.add(read);
        }

        final LocusIteratorByState locusIterator = new LocusIteratorByState(windowReads.iterator(), DownsamplingMethod.NONE, false, false, ReadUtils.getSamplesFromHeader(getHeaderForReads()), getHeaderForReads());
        final ActivityProfile activityProfile = new BandPassActivityProfile(null, 50, 0.002, BandPassActivityProfile.MAX_FILTER_SIZE, BandPassActivityProfile.DEFAULT_SIGMA, getHeaderForReads());

        List<AssemblyRegion> assemblyRegions = new ArrayList<>();
        for ( final AlignmentContext pileup : locusIterator ) {
            if ( ! activityProfile.isEmpty() && pileup.getLocation().getStart() != activityProfile.getEnd() + 1 ) {
                assemblyRegions.addAll(activityProfile.popReadyActiveRegions(100, 50, 300, true));
            }

            if ( window.getInterval().contains(pileup.getLocation()) ) {
                final SimpleInterval pileupInterval = new SimpleInterval(pileup.getLocation());
                final ReferenceBases refBase = new ReferenceBases(new byte[]{referenceContext.getBases()[pileup.getLocation().getStart() - referenceContext.getInterval().getStart()]}, pileupInterval);
                final ReferenceContext pileupRefContext = new ReferenceContext(new ReferenceMemorySource(refBase, getHeaderForReads().getSequenceDictionary()), pileupInterval);

                activityProfile.add(hcEngine.isActive(null, pileupRefContext, pileup));
            }
        }
        assemblyRegions.addAll(activityProfile.popReadyActiveRegions(100, 50, 300, true));

        final Iterator<GATKRead> reads = windowReads.iterator();
        AssemblyRegion lastRegion = null;
        for ( final AssemblyRegion region : assemblyRegions ) {
            if ( lastRegion != null ) {
                for ( final GATKRead lastRegionRead : lastRegion.getReads() ) {
                    if ( region.getExtendedSpan().overlaps(lastRegionRead) ) {
                        region.add(lastRegionRead);
                    }
                }
            }

            while ( reads.hasNext() ) {
                final GATKRead read = reads.next();
                if ( region.getExtendedSpan().overlaps(read) ) {
                    region.add(read);
                }
                else {
                    break;
                }
            }

            lastRegion = region;
        }

        return assemblyRegions;
    }

    @Override
    public Object onTraversalDone() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }

        if ( hcEngine != null ) {
            hcEngine.shutdown();
        }

        return null;
    }
}
