package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculator;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculators;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;


/**
 * Variant context utilities related to merging variant-context instances.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class ReferenceConfidenceVariantContextMerger {

    private static Comparable combineAnnotationValues( final List<Comparable> array ) {
        return MathUtils.median(array); // right now we take the median but other options could be explored
    }

    /**
     * Merges VariantContexts from gVCFs into a single hybrid.
     * Assumes that none of the input records are filtered.
     *
     * @param VCs     collection of unsorted genomic VCs
     * @param loc     the current location
     * @param refBase the reference allele to use if all contexts in the VC are spanning (i.e. don't start at the location in loc); if null, we'll return null in this case
     * @param removeNonRefSymbolicAllele if true, remove the <NON_REF> allele from the merged VC
     * @param samplesAreUniquified  if true, sample names have been uniquified
     * @return new VariantContext representing the merge of all VCs or null if it not relevant
     */
    public static VariantContext merge(final List<VariantContext> VCs, final GenomeLoc loc, final Byte refBase, final boolean removeNonRefSymbolicAllele,
                                       final boolean samplesAreUniquified) {
        // this can happen if e.g. you are using a dbSNP file that spans a region with no gVCFs
        if ( VCs == null || VCs.size() == 0 )
            return null;

        // establish the baseline info (sometimes from the first VC)
        final VariantContext first = VCs.get(0);
        final String name = first.getSource();

        // ref allele
        final Allele refAllele = determineReferenceAlleleGivenReferenceBase(VCs, loc, refBase);
        if ( refAllele == null )
            return null;

        // FinalAlleleSet contains the alleles of the new resulting VC
        // Using linked set in order to guarantee a stable order
        final LinkedHashSet<Allele> finalAlleleSet = new LinkedHashSet<>(10);
        // Reference goes first
        finalAlleleSet.add(refAllele);

        final Map<String, Object> attributes = new LinkedHashMap<>();
        final Set<String> rsIDs = new LinkedHashSet<>(1); // most of the time there's one id
        int depth = 0;
        final Map<String, List<Comparable>> annotationMap = new LinkedHashMap<>();
        final GenotypesContext genotypes = GenotypesContext.create();

        // In this list we hold the mapping of each variant context alleles.
        final List<Pair<VariantContext,List<Allele>>> vcAndNewAllelePairs = new ArrayList<>(VCs.size());
        // Keep track of whether we saw a spanning deletion and a non-spanning event
        boolean sawSpanningDeletion = false;
        boolean sawNonSpanningEvent = false;

        // cycle through and add info from the other VCs
        for ( final VariantContext vc : VCs ) {

            // if this context doesn't start at the current location then it must be a spanning event (deletion or ref block)
            final boolean isSpanningEvent = loc.getStart() != vc.getStart();
            // record whether it's also a spanning deletion/event (we know this because the VariantContext type is no
            // longer "symbolic" but "mixed" because there are real alleles mixed in with the symbolic non-ref allele)
            sawSpanningDeletion |= ( isSpanningEvent && vc.isMixed() ) || vc.getAlternateAlleles().contains(Allele.SPAN_DEL) ||
                    vc.getAlternateAlleles().contains(GATKVCFConstants.SPANNING_DELETION_SYMBOLIC_ALLELE_DEPRECATED );
            sawNonSpanningEvent |= ( !isSpanningEvent && vc.isMixed() );

            vcAndNewAllelePairs.add(new ImmutablePair<>(vc, isSpanningEvent ? replaceWithNoCallsAndDels(vc) : remapAlleles(vc, refAllele, finalAlleleSet)));
        }

        // Add <DEL> and <NON_REF> to the end if at all required in in the output.
        if ( sawSpanningDeletion && (sawNonSpanningEvent || !removeNonRefSymbolicAllele) ) finalAlleleSet.add(Allele.SPAN_DEL);
        if (!removeNonRefSymbolicAllele) finalAlleleSet.add(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);

        final List<Allele> allelesList = new ArrayList<>(finalAlleleSet);

        for ( final Pair<VariantContext,List<Allele>> pair : vcAndNewAllelePairs ) {
            final VariantContext vc = pair.getLeft();
            final List<Allele> remappedAlleles = pair.getRight();

            mergeRefConfidenceGenotypes(genotypes, vc, remappedAlleles, allelesList, samplesAreUniquified);

            // special case DP (add it up) for all events
            if ( vc.hasAttribute(VCFConstants.DEPTH_KEY) ) {
                depth += vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0);
            } else { // handle the gVCF case from the HaplotypeCaller
                for( final Genotype gt : vc.getGenotypes() ) {
                    depth += (gt.hasExtendedAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY) ? Integer.parseInt((String)gt.getAnyAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY)) : (gt.hasDP() ? gt.getDP() : 0));
                }
            }

            if ( loc.getStart() != vc.getStart() )
                continue;

            // special case ID (just preserve it)
            if ( vc.hasID() ) rsIDs.add(vc.getID());

            // add attributes
            addReferenceConfidenceAttributes(vc.getAttributes(), annotationMap);
        }

        // when combining annotations use the median value from all input VCs which had annotations provided
        for ( final Map.Entry<String, List<Comparable>> p : annotationMap.entrySet() ) {
            if ( ! p.getValue().isEmpty() ) {
                attributes.put(p.getKey(), combineAnnotationValues(p.getValue()));
            }
        }

        if ( depth > 0 )
            attributes.put(VCFConstants.DEPTH_KEY, String.valueOf(depth));

        // remove stale AC and AF based attributes
        removeStaleAttributesAfterMerge(attributes);

        final String ID = rsIDs.isEmpty() ? VCFConstants.EMPTY_ID_FIELD : Utils.join(",", rsIDs);

        // note that in order to calculate the end position, we need a list of alleles that doesn't include anything symbolic
        final VariantContextBuilder builder = new VariantContextBuilder().source(name).id(ID).alleles(allelesList)
                .chr(loc.getContig()).start(loc.getStart()).computeEndFromAlleles(nonSymbolicAlleles(allelesList), loc.getStart(), loc.getStart())
                .genotypes(genotypes).unfiltered().attributes(new TreeMap<>(attributes)).log10PError(CommonInfo.NO_LOG10_PERROR);  // we will need to re-genotype later

        return builder.make();
    }

    /**
     * @param list  the original alleles list
     * @return a non-null list of non-symbolic alleles
     */
    private static List<Allele> nonSymbolicAlleles(final List<Allele> list) {
        final List<Allele> result = new ArrayList<>(list.size());
        for ( final Allele allele : list ) {
            if ( !allele.isSymbolic() )
                result.add(allele);
        }
        return result;
    }

    /**
     * Determines the ref allele given the provided reference base at this position
     *
     * @param VCs     collection of unsorted genomic VCs
     * @param loc     the current location
     * @param refBase the reference allele to use if all contexts in the VC are spanning
     * @return new Allele or null if no reference allele/base is available
     */
    private static Allele determineReferenceAlleleGivenReferenceBase(final List<VariantContext> VCs, final GenomeLoc loc, final Byte refBase) {
        final Allele refAllele = GATKVariantContextUtils.determineReferenceAllele(VCs, loc);
        if ( refAllele == null )
            return ( refBase == null ? null : Allele.create(refBase, true) );
        return refAllele;
    }

    /**
     * Remove the stale attributes from the merged set
     *
     * @param attributes the attribute map
     */
    private static void removeStaleAttributesAfterMerge(final Map<String, Object> attributes) {
        attributes.remove(VCFConstants.ALLELE_COUNT_KEY);
        attributes.remove(VCFConstants.ALLELE_FREQUENCY_KEY);
        attributes.remove(VCFConstants.ALLELE_NUMBER_KEY);
        attributes.remove(GATKVCFConstants.MLE_ALLELE_COUNT_KEY);
        attributes.remove(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY);
        attributes.remove(VCFConstants.END_KEY);
    }

    /**
     * Adds attributes to the global map from the new context in a sophisticated manner
     *
     * @param myAttributes               attributes to add from
     * @param annotationMap              map of annotations for combining later
     */
    private static void addReferenceConfidenceAttributes(final Map<String, Object> myAttributes,
                                                         final Map<String, List<Comparable>> annotationMap) {
        for ( final Map.Entry<String, Object> p : myAttributes.entrySet() ) {
            final String key = p.getKey();
            final Object value = p.getValue();

            // add the annotation values to a list for combining later
            List<Comparable> values = annotationMap.get(key);
            if( values == null ) {
                values = new ArrayList<>();
                annotationMap.put(key, values);
            }
            try {
                final String stringValue = value.toString();
                // Branch to avoid unintentional, implicit type conversions that occur with the ? operator.
                if (stringValue.contains("."))
                    values.add(Double.parseDouble(stringValue));
                else
                    values.add(Integer.parseInt(stringValue));
            } catch (final NumberFormatException e) {
                // nothing to do
            }
        }
    }

    /**
     * This method does a couple of things:
     * <ul><li>
     *     remaps the vc alleles considering the differences between the final reference allele and its own reference,</li>
     * <li>
     *     collects alternative alleles present in variant context and add them to the {@code finalAlleles} set.
     * </li></ul>
     *
     * @param vc           the variant context.
     * @param refAllele    final reference allele.
     * @param finalAlleles where to add the final set of non-ref called alleles.
     * @return never {@code null}
     */
    //TODO as part of a larger refactoring effort {@link #remapAlleles} can be merged with {@link GATKVariantContextUtils#remapAlleles}.
    private static List<Allele> remapAlleles(final VariantContext vc, final Allele refAllele, final LinkedHashSet<Allele> finalAlleles) {

        final Allele vcRef = vc.getReference();
        final byte[] refBases = refAllele.getBases();
        final int extraBaseCount = refBases.length - vcRef.getBases().length;
        if (extraBaseCount < 0) throw new IllegalStateException("the wrong reference was selected");

        final List<Allele> result = new ArrayList<>(vc.getNAlleles());
        result.add(refAllele);

        for (final Allele a : vc.getAlternateAlleles()) {
            if (a.isSymbolic()) {
                result.add(a);
                // we always skip <NON_REF> when adding to finalAlleles; this is done outside if it applies.
                // we also skip <*:DEL> if there isn't a real alternate allele.
                if ( !a.equals(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE) && !vc.isSymbolic() )
                    finalAlleles.add(a);
            } else if ( a == Allele.SPAN_DEL ) {
                result.add(a);
                // we skip * if there isn't a real alternate allele.
                if ( !vc.isBiallelic() )
                    finalAlleles.add(a);
            } else if (a.isCalled()) {
                final Allele newAllele;
                if (extraBaseCount > 0) {
                    final byte[] oldBases = a.getBases();
                    final byte[] newBases = Arrays.copyOf(oldBases,oldBases.length + extraBaseCount);
                    System.arraycopy(refBases,refBases.length - extraBaseCount,newBases,oldBases.length,extraBaseCount);
                    newAllele = Allele.create(newBases,false);
                } else
                    newAllele = a;
                result.add(newAllele);
                finalAlleles.add(newAllele);
            } else { // NO_CALL and strange miscellanea
                result.add(a);
            }
        }
        return result;
    }

    /**
     * Replaces any alleles in the VariantContext with NO CALLS or the symbolic deletion allele as appropriate, except for the generic ALT allele
     *
     * @param vc   VariantContext with the alleles to replace
     * @return non-null list of alleles
     */
    private static List<Allele> replaceWithNoCallsAndDels(final VariantContext vc) {
        if ( vc == null ) throw new IllegalArgumentException("VariantContext cannot be null");

        final List<Allele> result = new ArrayList<>(vc.getNAlleles());

        // no-call the reference allele
        result.add(Allele.NO_CALL);

        // handle the alternate alleles
        for ( final Allele allele : vc.getAlternateAlleles() ) {
            final Allele replacement;
            if ( allele.equals(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE) )
                replacement = allele;
            else if ( allele.length() < vc.getReference().length() )
                replacement = Allele.SPAN_DEL;
            else
                replacement = Allele.NO_CALL;

            result.add(replacement);
        }
        return result;
    }

    /**
     * Merge into the context a new genotype represented by the given VariantContext for the provided list of target alleles.
     * This method assumes that none of the alleles in the VC overlaps with any of the alleles in the set.
     *
     * @param mergedGenotypes       the genotypes context to add to
     * @param VC                    the Variant Context for the sample
     * @param remappedAlleles       the list of remapped alleles for the sample
     * @param targetAlleles         the list of target alleles
     * @param samplesAreUniquified  true if sample names have been uniquified
     */
    private static void mergeRefConfidenceGenotypes(final GenotypesContext mergedGenotypes,
                                                    final VariantContext VC,
                                                    final List<Allele> remappedAlleles,
                                                    final List<Allele> targetAlleles,
                                                    final boolean samplesAreUniquified) {
        final int maximumPloidy = VC.getMaxPloidy(GATKVariantContextUtils.DEFAULT_PLOIDY);
        // the map is different depending on the ploidy, so in order to keep this method flexible (mixed ploidies)
        // we need to get a map done (lazily inside the loop) for each ploidy, up to the maximum possible.
        final int[][] genotypeIndexMapsByPloidy = new int[maximumPloidy + 1][];
        final int maximumAlleleCount = Math.max(remappedAlleles.size(),targetAlleles.size());
        int[] perSampleIndexesOfRelevantAlleles;

        for ( final Genotype g : VC.getGenotypes() ) {
            final String name;
            if (samplesAreUniquified)
                name = g.getSampleName() + "." + VC.getSource();
            else
                name = g.getSampleName();
            final int ploidy = g.getPloidy();
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(g).alleles(GATKVariantContextUtils.noCallAlleles(g.getPloidy()));
            genotypeBuilder.name(name);
            if (g.hasPL()) {
                // lazy initialization of the genotype index map by ploidy.
                perSampleIndexesOfRelevantAlleles = getIndexesOfRelevantAlleles(remappedAlleles, targetAlleles, VC.getStart(), g);
                final int[] genotypeIndexMapByPloidy = genotypeIndexMapsByPloidy[ploidy] == null
                        ? GenotypeLikelihoodCalculators.getInstance(ploidy, maximumAlleleCount).genotypeIndexMap(perSampleIndexesOfRelevantAlleles)
                        : genotypeIndexMapsByPloidy[ploidy];
                final int[] PLs = generatePL(g, genotypeIndexMapByPloidy);
                final int[] AD = g.hasAD() ? generateAD(g.getAD(), perSampleIndexesOfRelevantAlleles) : null;
                genotypeBuilder.PL(PLs).AD(AD);
            }
            mergedGenotypes.add(genotypeBuilder.make());
        }
    }

    /**
     * Composes a new likelihood array given the original genotype and the genotype index map.
     *
     * @param g the original genotype.
     * @param genotypeIndexMapByPloidy genotype index map. The ith element indicates what genotype in {@code g} corresponds
     *                                 to the ith genotype in the return likelihoods array.
     *
     * @throws NullPointerException if {@code g} or {@code genotypeIndexMapByPloidy} is {@code null}, or if {@code g}
     *    does not contain likelihoods.
     * @throws IndexOutOfBoundsException if {@code genotypeIndexMapByPloidy} contain non valid
     *  genotype indices given the likelihood array in {@code g}.
     *
     * @return never {@code null} but an array of exactly {@code genotypeIndexMapByPloidy.length} positions.
     */
    private static int[] generatePL(final Genotype g, final int[] genotypeIndexMapByPloidy) {
        final int[] PLs = new int[genotypeIndexMapByPloidy.length];
        final int[] oldPLs = g.getPL();
        for (int i = 0; i < PLs.length; i++)
            PLs[i] = oldPLs[genotypeIndexMapByPloidy[i]];
        return PLs;
    }

    /**
     * Determines the allele mapping from myAlleles to the targetAlleles, substituting the generic "<ALT>" as appropriate.
     * If the myAlleles set does not contain "<ALT>" as an allele, it throws an exception.
     *
     * @param remappedAlleles   the list of alleles to evaluate
     * @param targetAlleles     the target list of alleles
     * @param position          position to output error info
     * @param g                 genotype from which targetAlleles are derived
     * @return non-null array of ints representing indexes
     */
    protected static int[] getIndexesOfRelevantAlleles(final List<Allele> remappedAlleles, final List<Allele> targetAlleles, final int position, final Genotype g) {

        if ( remappedAlleles == null || remappedAlleles.size() == 0 ) throw new IllegalArgumentException("The list of input alleles must not be null or empty");
        if ( targetAlleles == null || targetAlleles.size() == 0 ) throw new IllegalArgumentException("The list of target alleles must not be null or empty");

        if ( !remappedAlleles.contains(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE) )
            throw new UserException("The list of input alleles must contain " + GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE + " as an allele but that is not the case at position " + position + "; please use the Haplotype Caller with gVCF output to generate appropriate records");

        final int indexOfNonRef = remappedAlleles.indexOf(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
        final int[] indexMapping = new int[targetAlleles.size()];

        // the reference likelihoods should always map to each other (even if the alleles don't)
        indexMapping[0] = 0;

        // create the index mapping, using the <NON-REF> allele whenever such a mapping doesn't exist
        for ( int i = 1; i < targetAlleles.size(); i++ ) {
            final Allele targetAllele = targetAlleles.get(i);

            // if there’s more than 1 DEL allele then we need to use the best one
            if ( targetAllele == Allele.SPAN_DEL && g.hasPL() ) {
                final int occurrences = Collections.frequency(remappedAlleles, Allele.SPAN_DEL);
                if ( occurrences > 1 ) {
                    final int indexOfBestDel = indexOfBestDel(remappedAlleles, g.getPL(), g.getPloidy());
                    indexMapping[i] = ( indexOfBestDel == -1 ? indexOfNonRef : indexOfBestDel );
                    continue;
                }
            }

            final int indexOfRemappedAllele = remappedAlleles.indexOf(targetAllele);
            indexMapping[i] = indexOfRemappedAllele == -1 ? indexOfNonRef : indexOfRemappedAllele;
        }

        return indexMapping;
    }

    /**
     * Returns the index of the best spanning deletion allele based on AD counts
     *
     * @param alleles   the list of alleles
     * @param PLs       the list of corresponding PL values
     * @param ploidy    the ploidy of the sample
     * @return the best index or -1 if not found
     */
    private static int indexOfBestDel(final List<Allele> alleles, final int[] PLs, final int ploidy) {
        int bestIndex = -1;
        int bestPL = Integer.MAX_VALUE;

        for ( int i = 0; i < alleles.size(); i++ ) {
            if ( alleles.get(i) == Allele.SPAN_DEL ) {
                final int homAltIndex = findHomIndex(i, ploidy, alleles.size());
                final int PL = PLs[homAltIndex];
                if ( PL < bestPL ) {
                    bestIndex = i;
                    bestPL = PL;
                }
            }
        }

        return bestIndex;
    }

    /**
     * Returns the index of the PL that represents the homozygous genotype of the given i'th allele
     *
     * @param i           the index of the allele with the list of alleles
     * @param ploidy      the ploidy of the sample
     * @param numAlleles  the total number of alleles
     * @return the hom index
     */
    private static int findHomIndex(final int i, final int ploidy, final int numAlleles) {
        // some quick optimizations for the common case
        if ( ploidy == 2 )
            return i + (i * (i + 1) / 2);   // this is straight from the VCF spec on PLs
        if ( ploidy == 1 )
            return i;

        final GenotypeLikelihoodCalculator calculator = GenotypeLikelihoodCalculators.getInstance(ploidy, numAlleles);
        final int[] alleleIndexes = new int[ploidy];
        Arrays.fill(alleleIndexes, i);
        return calculator.allelesToIndex(alleleIndexes);
    }

    /**
     * Generates a new AD array by adding zeros for missing alleles given the set of indexes of the Genotype's current
     * alleles from the original AD.
     *
     * @param originalAD    the original AD to extend
     * @param indexesOfRelevantAlleles the indexes of the original alleles corresponding to the new alleles
     * @return non-null array of new AD values
     */
    protected static int[] generateAD(final int[] originalAD, final int[] indexesOfRelevantAlleles) {
        if ( originalAD == null || indexesOfRelevantAlleles == null ) throw new IllegalArgumentException("The list of input AD values and alleles must not be null");

        final int numADs = indexesOfRelevantAlleles.length;
        final int[] newAD = new int[numADs];

        for ( int i = 0; i < numADs; i++ ) {
            final int oldIndex = indexesOfRelevantAlleles[i];
            if ( oldIndex >= originalAD.length )
                newAD[i] = 0;
            else
                newAD[i] = originalAD[oldIndex];
        }

        return newAD;
    }
}
