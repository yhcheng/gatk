package org.broadinstitute.hellbender.tools.walkers.variantutils;

//------------------------------------------------------
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.tribble.Feature;
//import htsjdk.samtools.util.Locatable;

import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.utils.SampleUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.commandline.*;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
//import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;

import java.io.File;
import java.util.*;

/**
 * Combine per-sample gVCF files produced by HaplotypeCaller into a multi-sample gVCF file
 *
 * <p>
 * CombineGVCFs is meant to be used for hierarchical merging of gVCFs that will eventually be input into GenotypeGVCFs.
 * One would use this tool when needing to genotype too large a number of individual gVCFs; instead of passing them
 * all in to GenotypeGVCFs, one would first use CombineGVCFs on smaller batches of samples and then pass these combined
 * gVCFs to GenotypeGVCFs.</p>
 *
 * <h3>Input</h3>
 * <p>
 * Two or more Haplotype Caller gVCFs to combine.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A combined multisample gVCF.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T CombineGVCFs \
 *   -R reference.fasta \
 *   --variant sample1.g.vcf \
 *   --variant sample2.g.vcf \
 *   -o cohort.g.vcf
 * </pre>
 *
 * <h3>Caveat</h3>
 * <p>Only gVCF files produced by HaplotypeCaller (or CombineGVCFs) can be used as input for this tool. Some other
 * programs produce files that they call gVCFs but those lack some important information (accurate genotype likelihoods
 * for every position) that GenotypeGVCFs requires for its operation.</p>
 *
 */
//@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARMANIP, extraDocs = {CommandLineGATK.class} )
//@Reference(window=@Window(start=0,stop=1))
//public class CombineGVCFs extends RodWalker<CombineGVCFs.PositionalState, CombineGVCFs.OverallState> {
@CommandLineProgramProperties(
                summary = "Combine GVCF files produced by HaplotypCaller into multi-sample GVCF file",
                oneLineSummary = "Combine GVCF files",
                programGroup = VariantProgramGroup.class
)

public class CombineGVCFs extends MultiVariantWalker {

    protected final class PositionalState {
        final List<VariantContext> VCs;
        final Set<String> samples = new HashSet<>();
        final byte[] refBases;
        //final GenomeLoc loc;
        final GenomeLoc loc;
        public PositionalState(final List<VariantContext> VCs, final byte[] refBases, final GenomeLoc loc) {
            this.VCs = VCs;
            for(final VariantContext vc : VCs){
                samples.addAll(vc.getSampleNames());
            }
            this.refBases = refBases;
            this.loc = loc;
        }
    }

    protected final class OverallState {
        final LinkedList<VariantContext> VCs = new LinkedList<>();
        final Set<String> samples = new HashSet<>();
        GenomeLoc prevPos = null;
        byte refAfterPrevPos;

        public OverallState() {}
    }

    //final private List<VariantContext> variants = new ArrayList<>();

    @Argument(doc="File to which the combined gVCF should be written")
    private File outFile = null;

    protected VariantContextWriter vcfWriter = null;

    @Argument(fullName="convertToBasePairResolution", shortName="bpResolution", doc = "If specified, convert banded gVCFs to all-sites gVCFs", optional=true)
    protected boolean USE_BP_RESOLUTION = false;

    /**
     * To reduce file sizes our gVCFs group similar reference positions into bands.  However, there are cases when users will want to know that no bands
     * span across a given genomic position (e.g. when scatter-gathering jobs across a compute farm).  The option below enables users to break bands at
     * pre-defined positions.  For example, a value of 10,000 would mean that we would ensure that no bands span across chr1:10000, chr1:20000, etc.
     *
     * Note that the --convertToBasePairResolution argument is just a special case of this argument with a value of 1.
     */
    @Argument(fullName="breakBandsAtMultiplesOf", shortName="breakBandsAtMultiplesOf", doc = "If > 0, reference bands will be broken up at genomic positions that are multiples of this number", optional=true)
    protected int multipleAtWhichToBreakBands = 0;

    private GenomeLocParser genomeLocParser;

    //private PositionalState posState;
    private OverallState overallState = new OverallState();

    @Override
    public void onTraversalStart() {
        final List<FeatureInput<VariantContext>> featureInputs = getDrivingVariantsFeatureInputs();
        //FeatureInput<VariantContext> vc = variantCollections.get(0);

        // take care of the VCF headers
        final Map<String, VCFHeader> vcfRods = new HashMap<String, VCFHeader>();
        List<VCFHeader> vcfHeaders= getHeadersForVariants();
        for (int i = 0; i < featureInputs.size(); i++) {
            FeatureInput<VariantContext> featureInput = featureInputs.get(i);
            vcfRods.put(featureInput.getName(), vcfHeaders.get(i));
        }

        // Initialize VCF header lines
        final Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), true);
        // TODO - should this header be added ? GATK doesn't add one for CombineGVCfs but does for SelectVariants
        //headerLines.add(new VCFHeaderLine("source", "CombineGVCFs"));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));   // needed for gVCFs without DP tags

        final Set<String> samples = SampleUtils.getSampleList(vcfRods, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);

        final File refFile = referenceArguments.getReferenceFile();
        SAMSequenceDictionary samDictionary = this.getBestAvailableSequenceDictionary();
        Set<VCFHeaderLine> actualLines = null;
        if (refFile != null && samDictionary!= null) {
            // TOOD: suppressReferencelines ?
            //actualLines = withUpdatedContigsAsLines(headerLines, refFile, samDictionary, suppressReferenceLines);
            actualLines = withUpdatedContigsAsLines(headerLines, refFile, samDictionary, false);
        }
        else {
            actualLines = headerLines;
        }

        VariantContextWriterBuilder vcWriterBuilder =
                new VariantContextWriterBuilder()
                    .setOutputFile(outFile)
                    .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
                     //.unsetOption(Options.WRITE_FULL_FORMAT_FIELD)
                    .unsetOption(Options.INDEX_ON_THE_FLY);
        //if (lenientVCFProcesssing)
        //  vcWriterBuilder = vcWriterBuilder.setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER);
        vcfWriter = vcWriterBuilder.build();
        vcfWriter.writeHeader(new VCFHeader(actualLines, samples));

        // collect the actual rod bindings into a list for use later
        //for ( final RodBindingCollection<VariantContext> variantCollection : variantCollections )
        //    variants.addAll(variantCollection.getRodBindings());
        //variants.add()
        //for ( final FeatureInput<VariantContext> fi : variantCollections )
        //    variants.add(fi.);

        ReferenceSequenceFile rsf = ReferenceSequenceFileFactory.getReferenceSequenceFile(referenceArguments.getReferenceFile());
        genomeLocParser = new GenomeLocParser(rsf);

        // optimization to prevent mods when we always just want to break bands
        if ( multipleAtWhichToBreakBands == 1 )
            USE_BP_RESOLUTION = true;
    }

    @Override
    public void apply(VariantContext vc, ReadsContext readsContext, ReferenceContext ref, FeatureContext featureContext) {
        // TODO - is this equvalent to getReferenceIndex ?
        int contigIndex = readsContext.getInterval().getStart();

        final GenomeLoc initLoc = genomeLocParser.createGenomeLoc(vc);
        int leadSize = initLoc.getStop() - initLoc.getStart();
        //ref.setWindow(0, leadSize);
        byte bases[] = ref.getBases();
        for (int i = 0; i < leadSize + 1; i ++) {
            // TODO is this right ? vc.getStart ?
            //return new PositionalState(tracker.getValues(variants, loc), ref.getBases(), loc);
            // there MUST be a better way..
            GenomeLoc loc = initLoc.setStart(initLoc, initLoc.getStart() + i);
            loc = initLoc.setStop(loc, initLoc.getStart() + i);
            if (i != 0) {
                int len = bases.length;
                bases = Arrays.copyOfRange(bases, 1, bases.length);
            }
            PositionalState posState = new PositionalState(featureContext.getValues(getDrivingVariantsFeatureInputs(), loc.getStart()), bases, loc);

            // startingState == posState
            // previousState == overallState
            // overallState = newState
            if ( !posState.VCs.isEmpty() ) {
                if ( ! okayToSkipThisSite(posState, overallState) ) {
                    // TODO: this is way ugly
                    GenomeLoc newLoc = posState.loc.setStart(posState.loc, posState.loc.getStart() - 1);
                    newLoc = newLoc.setStop(newLoc, posState.loc.getStop() - 1);
                    endPreviousStates(overallState, newLoc, posState, false);
                }
                overallState.VCs.addAll(posState.VCs);
                for(final VariantContext vContext : overallState.VCs){
                    overallState.samples.addAll(vContext.getSampleNames());
                }
            }

            if ( breakBand(posState.loc) || containsEndingContext(overallState.VCs, posState.loc.getStart()) ) {
                endPreviousStates(overallState, posState.loc, posState, true);
            }
        }
    }

    @Override
    public Object onTraversalDone() {
        try {
            return null;
        } finally {
            vcfWriter.close();
        }
    }

    private static Set<VCFHeaderLine> withUpdatedContigsAsLines(
            final Set<VCFHeaderLine> oldLines,
            final File referenceFile,
            final SAMSequenceDictionary refDict,
            final boolean referenceNameOnly) {
        final Set<VCFHeaderLine> lines = new LinkedHashSet<VCFHeaderLine>(oldLines.size());

        for (final VCFHeaderLine line : oldLines) {
            if (line instanceof VCFContigHeaderLine)
                continue; // skip old contig lines
            if (line.getKey().equals(VCFHeader.REFERENCE_KEY))
                continue; // skip the old reference key
            lines.add(line);
        }

        for (final VCFHeaderLine contigLine : makeContigHeaderLines(refDict, referenceFile))
            lines.add(contigLine);

        final String referenceValue;
        if (referenceFile != null) {
            if (referenceNameOnly) {
                final int extensionStart = referenceFile.getName().lastIndexOf(".");
                referenceValue = extensionStart == -1 ? referenceFile.getName() : referenceFile.getName().substring(0, extensionStart);
            }
            else {
                referenceValue = "file://" + referenceFile.getAbsolutePath();
            }
            lines.add(new VCFHeaderLine(VCFHeader.REFERENCE_KEY, referenceValue));
        }
        return lines;
    }

    /**
     * Create VCFHeaderLines for each refDict entry, and optionally the assembly if referenceFile != null
     * @param refDict reference dictionary
     * @param referenceFile for assembly name.  May be null
     * @return list of vcf contig header lines
     */
    public static List<VCFContigHeaderLine> makeContigHeaderLines(final SAMSequenceDictionary refDict,
                                                                  final File referenceFile) {
        final List<VCFContigHeaderLine> lines = new ArrayList<VCFContigHeaderLine>();
        final String assembly = referenceFile != null ? getReferenceAssembly(referenceFile.getName()) : null;
        for (final SAMSequenceRecord contig : refDict.getSequences())
            lines.add(makeContigHeaderLine(contig, assembly));
        return lines;
    }

    private static VCFContigHeaderLine makeContigHeaderLine(final SAMSequenceRecord contig, final String assembly) {
        final Map<String, String> map = new LinkedHashMap<String, String>(3);
        map.put("ID", contig.getSequenceName());
        map.put("length", String.valueOf(contig.getSequenceLength()));
        if (assembly != null) {
            map.put("assembly", assembly);
        }
        return new VCFContigHeaderLine(map, contig.getSequenceIndex());
    }

    // TODO: the GATK code does this to include in the headers - seems pretty sketchy...
    private static String getReferenceAssembly(final String refPath) {
        // This doesn't need to be perfect as it's not a required VCF header line, but we might as well give it a shot
        String assembly = null;
        if (refPath.contains("b37") || refPath.contains("v37"))
            assembly = "b37";
        else if (refPath.contains("b36"))
            assembly = "b36";
        else if (refPath.contains("hg18"))
            assembly = "hg18";
        else if (refPath.contains("hg19"))
            assembly = "hg19";
        return assembly;
    }

    /*
    public void onTraversalDone(final OverallState state) {
        // there shouldn't be any state left unless the user cut in the middle of a gVCF block
        if ( !state.VCs.isEmpty() )
            logger.warn("You have asked for an interval that cuts in the middle of one or more gVCF blocks. Please note that this will cause you to lose records that don't end within your interval.");
    }

    public PositionalState map(final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context) {
        if ( tracker == null ) // RodWalkers can make funky map calls
            return null;

        final GenomeLoc loc = ref.getLocus();
        return new PositionalState(tracker.getValues(variants, loc), ref.getBases(), loc);
    }

    public OverallState reduceInit() {
        return new OverallState();
    }

    public OverallState reduce(final PositionalState startingStates, final OverallState previousState) {
        if ( startingStates == null )
            return previousState;

        if ( !startingStates.VCs.isEmpty() ) {
            if ( ! okayToSkipThisSite(startingStates, previousState) )
                endPreviousStates(previousState, startingStates.loc.incPos(-1), startingStates, false);
            previousState.VCs.addAll(startingStates.VCs);
            for(final VariantContext vc : previousState.VCs){
                previousState.samples.addAll(vc.getSampleNames());
            }

        }

        if ( breakBand(startingStates.loc) || containsEndingContext(previousState.VCs, startingStates.loc.getStart()) ) {
            endPreviousStates(previousState, startingStates.loc, startingStates, true);
        }

        return previousState;
    }
    */

    /**
     * Should we break bands at the given position?
     *
     * @param loc  the genomic location to evaluate against
     *
     * @return true if we should ensure that bands should be broken at the given position, false otherwise
     */
    private boolean breakBand(final GenomeLoc loc) {
        return USE_BP_RESOLUTION ||
                (loc != null && multipleAtWhichToBreakBands > 0 && (loc.getStart()+1) % multipleAtWhichToBreakBands == 0);  // add +1 to the loc because we want to break BEFORE this base
    }

    /**
     * Is it okay to skip the given position?
     *
     * @param startingStates  state information for this position
     * @param previousState   state information for the last position for which we created a VariantContext
     * @return true if it is okay to skip this position, false otherwise
     */
    private boolean okayToSkipThisSite(final PositionalState startingStates, final OverallState previousState) {
        final int thisPos = startingStates.loc.getStart();
        final GenomeLoc lastPosRun = previousState.prevPos;
        Set<String> intersection = new HashSet<String>(startingStates.samples);
        intersection.retainAll(previousState.samples);

        //if there's a starting VC with a sample that's already in a current VC, don't skip this position
        return lastPosRun != null && thisPos == lastPosRun.getStart() + 1 && intersection.isEmpty();
    }

    /**
     * Does the given list of VariantContexts contain any whose context ends at the given position?
     *
     * @param VCs  list of VariantContexts
     * @param pos  the position to check against
     * @return true if there are one or more VCs that end at pos, false otherwise
     */
    private boolean containsEndingContext(final List<VariantContext> VCs, final int pos) {
        if ( VCs == null ) throw new IllegalArgumentException("The list of VariantContexts cannot be null");

        for ( final VariantContext vc : VCs ) {
            if ( isEndingContext(vc, pos) )
                return true;
        }
        return false;
    }

    /**
     * Does the given variant context end (in terms of reference blocks, not necessarily formally) at the given position.
     * Note that for the purposes of this method/tool, deletions are considered to be single base events (as opposed to
     * reference blocks), hence the check for the number of alleles (because we know there will always be a <NON_REF> allele).
     *
     * @param vc   the variant context
     * @param pos  the position to query against
     * @return true if this variant context "ends" at this position, false otherwise
     */
    private boolean isEndingContext(final VariantContext vc, final int pos) {
        return vc.getNAlleles() > 2 || vc.getEnd() == pos;
    }

    /**
     * Disrupt the VariantContexts so that they all stop at the given pos, write them out, and put the remainder back in the list.
     * @param state   the previous state with list of active VariantContexts
     * @param pos   the position for the starting VCs
     * @param startingStates the state for the starting VCs
     * @param atCurrentPosition  indicates whether we output a variant at the current position, independent of VCF start/end, i.e. in BP resolution mode
     */
    private void endPreviousStates(final OverallState state, final GenomeLoc pos, final PositionalState startingStates, boolean atCurrentPosition) {

        final byte refBase = startingStates.refBases[0];
        //if we're in BP resolution mode or a VC ends at the current position then the reference for the next output VC (refNextBase)
        // will be advanced one base
        final byte refNextBase = (atCurrentPosition) ? (startingStates.refBases.length > 1 ? startingStates.refBases[1] : (byte)'N' ): refBase;

        final List<VariantContext> stoppedVCs = new ArrayList<>(state.VCs.size());

        for ( int i = state.VCs.size() - 1; i >= 0; i-- ) {
            final VariantContext vc = state.VCs.get(i);
            //the VC for the previous state will be stopped if its position is previous to the current position or it we've moved to a new contig
        // TODO is getChr -> getContig right ?
            //if ( vc.getStart() <= pos.getStart() || !vc.getChr().equals(pos.getContig())) {
            if ( vc.getStart() <= pos.getStart() || !vc.getContig().equals(pos.getContig())) {

                stoppedVCs.add(vc);

                // if it was ending anyways, then remove it from the future state
                if ( vc.getEnd() == pos.getStart()) {
                    state.samples.removeAll(vc.getSampleNames());
                    state.VCs.remove(i);
                    continue; //don't try to remove twice
                }

                //if ending vc is the same sample as a starting VC, then remove it from the future state
                if(startingStates.VCs.size() > 0 && !atCurrentPosition && startingStates.samples.containsAll(vc.getSampleNames())) {
                    state.samples.removeAll(vc.getSampleNames());
                    state.VCs.remove(i);
                }
            }
        }

        //output the stopped VCs if there is no previous output (state.prevPos == null) or our current position is past
        // the last write position (state.prevPos)
        //NOTE: BP resolution with have current position == state.prevPos because it gets output via a different control flow
        if ( !stoppedVCs.isEmpty() &&  (state.prevPos == null || pos.isPast(state.prevPos) )) {
            // TODO is getChr -> getContig right ?
            //final GenomeLoc gLoc = genomeLocParser.createGenomeLoc(stoppedVCs.get(0).getChr(), pos.getStart());
            final GenomeLoc gLoc = genomeLocParser.createGenomeLoc(stoppedVCs.get(0).getContig(), pos.getStart());

            // we need the specialized merge if the site contains anything other than ref blocks
            final VariantContext mergedVC;
            if ( containsTrueAltAllele(stoppedVCs) ) {
                // TODO
                mergedVC = ReferenceConfidenceVariantContextMerger.merge(stoppedVCs, gLoc, refBase, false, false);
            } else {
                mergedVC = referenceBlockMerge(stoppedVCs, state, pos.getStart());
            }

            vcfWriter.add(mergedVC);
            state.prevPos = gLoc;
            state.refAfterPrevPos = refNextBase;
        }
    }

    /**
     * Combine a list of reference block VariantContexts.
     * We can't use GATKVariantContextUtils.simpleMerge() because it is just too slow for this sort of thing.
     *
     * @param VCs   the variant contexts to merge
     * @param state the state object
     * @param end   the end of this block (inclusive)
     * @return a new merged VariantContext
     */
    private VariantContext referenceBlockMerge(final List<VariantContext> VCs, final OverallState state, final int end) {

        final VariantContext first = VCs.get(0);

        // ref allele and start
        final Allele refAllele;
        final int start;
    // TODO is getChr->getContig right >
        //if ( state.prevPos == null || !state.prevPos.getContig().equals(first.getChr()) || first.getStart() >= state.prevPos.getStart() + 1) {
        if ( state.prevPos == null || !state.prevPos.getContig().equals(first.getContig()) || first.getStart() >= state.prevPos.getStart() + 1) {
            start = first.getStart();
            refAllele = first.getReference();
        } else {
            start = state.prevPos.getStart() + 1;
            refAllele = Allele.create(state.refAfterPrevPos, true);
        }

        // attributes
        final Map<String, Object> attrs = new HashMap<>(1);
        if ( !USE_BP_RESOLUTION && end != start )
            attrs.put(VCFConstants.END_KEY, Integer.toString(end));

        // genotypes
        final GenotypesContext genotypes = GenotypesContext.create();
        for ( final VariantContext vc : VCs ) {
            for ( final Genotype g : vc.getGenotypes() )
                genotypes.add(new GenotypeBuilder(g).alleles(GATKVariantContextUtils.noCallAlleles(g.getPloidy())).make());
        }

        // TODO: move to GATKVCFConstants
        final String NON_REF_SYMBOLIC_ALLELE_NAME = "NON_REF";
        final Allele NON_REF_SYMBOLIC_ALLELE = Allele.create("<"+NON_REF_SYMBOLIC_ALLELE_NAME+">", false); // represents any possible non-ref allele at this site

    //  TODO is get Chr -> getContig right ?
        //return new VariantContextBuilder("", first.getChr(), start, end, Arrays.asList(refAllele, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE)).attributes(attrs).genotypes(genotypes).make();
        //return new VariantContextBuilder("", first.getChr(), start, end, Arrays.asList(refAllele, NON_REF_SYMBOLIC_ALLELE)).attributes(attrs).genotypes(genotypes).make();
        return new VariantContextBuilder("", first.getContig(), start, end, Arrays.asList(refAllele, NON_REF_SYMBOLIC_ALLELE)).attributes(attrs).genotypes(genotypes).make();
    }

    /**
     * Does the given list of VariantContexts contain any with an alternate allele other than <NON_REF>?
     *
     * @param VCs  list of VariantContexts
     * @return true if there are one or more VCs that contain a true alternate allele, false otherwise
     */
    private boolean containsTrueAltAllele(final List<VariantContext> VCs) {
        if ( VCs == null ) throw new IllegalArgumentException("The list of VariantContexts cannot be null");

        for ( final VariantContext vc : VCs ) {
            if ( vc.getNAlleles() > 2 )
                return true;
        }
        return false;
    }
}
