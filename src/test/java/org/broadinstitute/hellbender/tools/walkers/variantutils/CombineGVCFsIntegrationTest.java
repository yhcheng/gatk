package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.VariantContext;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;


public class CombineGVCFsIntegrationTest extends CommandLineProgramTest {

    public static final String b37KGReference = "/humgen/1kg/reference/human_g1k_v37.fasta";

    public String baseTestString(String args) {
        return " -R " + b37KGReference + " -o %s -V "
                + getToolTestDataDir() + "gvcfExample1.vcf -V " + getToolTestDataDir() + "gvcfExample2.vcf" + " " + args;
    }

    @Test
    public void testOneStartsBeforeTwoAndEndsAfterwards() throws Exception {
        final String expectedFile = getToolTestDataDir() + "expected/" + "OneStartsBeforeTwoAndEndsAfterwards.vcf";
        final String cmd = baseTestString(" -L 1:69485-69509");

        final IntegrationTestSpec spec = new IntegrationTestSpec (cmd, Arrays.asList(expectedFile));

        spec.executeTest("testOneStartsBeforeTwoAndEndsAfterwards", this);

        final File gVCF = new File(expectedFile);
        final List<VariantContext> allVCs = GATKVariantContextUtils.readVCF(gVCF).getRight();

        Assert.assertEquals(allVCs.size(), 2, "Observed: " + allVCs);

        final VariantContext first = allVCs.get(0);
        Assert.assertEquals(first.getStart(), 69491);
        Assert.assertEquals(first.getEnd(), 69497);
        Assert.assertEquals(first.getGenotypes().size(), 2);
        Assert.assertTrue(first.getGenotype("NA1").isNoCall());
        Assert.assertTrue(first.getGenotype("NA2").isNoCall());

        final VariantContext second = allVCs.get(1);
        Assert.assertEquals(second.getStart(), 69498);
        Assert.assertEquals(second.getEnd(), 69506);
        Assert.assertEquals(second.getGenotypes().size(), 2);
        Assert.assertTrue(second.getGenotype("NA1").isNoCall());
        Assert.assertTrue(second.getGenotype("NA2").isNoCall());
    }

    @Test(enabled = true)
    public void testTetraploidRun() throws IOException {
        final String expectedFile = getToolTestDataDir() + "expected/" + "combineSingleSamplePipelineGVCFTest.vcf";
        IntegrationTestSpec  spec = new IntegrationTestSpec (
                " -R " + b37KGReference + " -o %s " + "-V:sample1 " + getToolTestDataDir() + "tetraploid-gvcf-1.vcf" +
                        " -V:sample2 " + getToolTestDataDir() + "tetraploid-gvcf-2.vcf" +
                        " -V:sample3 " + getToolTestDataDir() + "tetraploid-gvcf-3.vcf" +
                        " -L " + getToolTestDataDir() + "tetraploid-gvcfs.intervals",
                Arrays.asList(expectedFile));
                //Arrays.asList("7b3153135e4f8e1d137d3f4beb46f182"));
        spec.executeTest("combineSingleSamplePipelineGVCFTest", this);
    }

    @Test(enabled= true)
    public void testMixedPloidyRun() throws IOException {
        final String expectedFile = getToolTestDataDir() + "expected/" + "combineSingleSamplePipelineGVCFMixed.vcf";
        IntegrationTestSpec  spec = new IntegrationTestSpec (
                " -R " + b37KGReference + " -o %s -V:sample1 " + getToolTestDataDir() + "haploid-gvcf-1.vcf" +
                        " -V:sample2 " + getToolTestDataDir() + "tetraploid-gvcf-2.vcf" +
                        " -V:sample3 " + getToolTestDataDir() + "diploid-gvcf-3.vcf" +
                        " -L " + getToolTestDataDir() + "tetraploid-gvcfs.intervals",
                Arrays.asList(expectedFile));
                //Arrays.asList("4f546634213ece6f08ec9258620b92bb"));
        spec.executeTest("combineSingleSamplePipelineGVCFMixed", this);
    }

    @Test
    public void testTwoSpansManyBlocksInOne() throws Exception {
        final String expectedFile = getToolTestDataDir() + "expected/" + "TwoSpansManyBlocksInOne.vcf";
        final String cmd = baseTestString(" -L 1:69512-69634");
        final IntegrationTestSpec spec = new IntegrationTestSpec (cmd, Arrays.asList(expectedFile));
        spec.executeTest("TwoSpansManyBlocksInOne", this);
        final File gVCF = new File(expectedFile);
        final List<VariantContext> allVCs = GATKVariantContextUtils.readVCF(gVCF).getRight();

        Assert.assertEquals(allVCs.size(), 5);
    }

    @Test
    public void testOneHasAltAndTwoHasNothing() throws Exception {
        final String expectedFile = getToolTestDataDir() + "expected/" + "OneHasAltAndTwoHasNothing.vcf";
        final String cmd = baseTestString(" -L 1:69511");
        final IntegrationTestSpec spec = new IntegrationTestSpec (cmd, Arrays.asList(expectedFile));
        spec.executeTest("OneHasAltAndTwoHasNothing", this);
        final File gVCF = new File(expectedFile);
        final List<VariantContext> allVCs = GATKVariantContextUtils.readVCF(gVCF).getRight();

        Assert.assertEquals(allVCs.size(), 1);

        final VariantContext first = allVCs.get(0);
        Assert.assertEquals(first.getStart(), 69511);
        Assert.assertEquals(first.getEnd(), 69511);
        Assert.assertEquals(first.getGenotypes().size(), 2);
    }

    @Test
    public void testOneHasAltAndTwoHasRefBlock() throws Exception {
        final String expectedFile = getToolTestDataDir() + "expected/" + "OneHasAltAndTwoHasRefBlock.vcf";
        final String cmd = baseTestString(" -L 1:69635");
        final IntegrationTestSpec  spec = new IntegrationTestSpec (cmd, Arrays.asList(expectedFile));
        spec.executeTest("testOneHasAltAndTwoHasRefBlock", this);
        final File gVCF = new File(expectedFile);
        final List<VariantContext> allVCs = GATKVariantContextUtils.readVCF(gVCF).getRight();

        Assert.assertEquals(allVCs.size(), 1);

        final VariantContext first = allVCs.get(0);
        Assert.assertEquals(first.getStart(), 69635);
        Assert.assertEquals(first.getEnd(), 69635);
        Assert.assertEquals(first.getNAlleles(), 3);
        Assert.assertEquals(first.getGenotypes().size(), 2);
    }

    @Test
    public void testOneHasDeletionAndTwoHasRefBlock() throws Exception {
        final String expectedFile = getToolTestDataDir() + "expected/" + "OneHasDeletionAndTwoHasRefBlock.vcf";
        final String cmd = baseTestString(" -L 1:69772-69783");
        final IntegrationTestSpec  spec = new IntegrationTestSpec (cmd, Arrays.asList(expectedFile));
        spec.executeTest("testOneHasDeletionAndTwoHasRefBlock", this);
        final File gVCF = new File(expectedFile);
        final List<VariantContext> allVCs = GATKVariantContextUtils.readVCF(gVCF).getRight();

        Assert.assertEquals(allVCs.size(), 3);

        final VariantContext first = allVCs.get(0);
        Assert.assertEquals(first.getStart(), 69772);
        Assert.assertEquals(first.getEnd(), 69776);
        Assert.assertEquals(first.getNAlleles(), 3);
        Assert.assertEquals(first.getGenotypes().size(), 2);

        final VariantContext second = allVCs.get(1);
        Assert.assertEquals(second.getStart(), 69773);
        Assert.assertEquals(second.getEnd(), 69774);
        Assert.assertEquals(second.getGenotypes().size(), 2);

        final VariantContext third = allVCs.get(2);
        Assert.assertEquals(third.getStart(), 69775);
        Assert.assertEquals(third.getEnd(), 69783);
        Assert.assertEquals(third.getGenotypes().size(), 2);
    }

    @Test
    public void testMD5s() throws Exception {
        final String expectedFile = getToolTestDataDir() + "expected/" + "testMD5s.vcf";
        final String cmd = baseTestString(" -L 1:69485-69791");
        final IntegrationTestSpec  spec = new IntegrationTestSpec (cmd, Arrays.asList(expectedFile));
        spec.executeTest("testMD5s", this);
    }

    @Test
    public void testBasepairResolutionOutput() throws Exception {
        final String expectedFile = getToolTestDataDir() + "expected/" + "BasepairResolutionOutput.vcf";
        final String cmd = baseTestString(" -L 1:69485-69791 --convertToBasePairResolution");
        final IntegrationTestSpec  spec = new IntegrationTestSpec (cmd, Arrays.asList(expectedFile));
        spec.executeTest("testBasepairResolutionOutput", this);
    }

    @Test
    public void testBreakBlocks() throws Exception {
        final String expectedFile = getToolTestDataDir() + "expected/" + "BreakBlocks.vcf";
        final String cmd = baseTestString(" -L 1:69485-69791 --breakBandsAtMultiplesOf 5");
        final IntegrationTestSpec  spec = new IntegrationTestSpec (cmd, Arrays.asList(expectedFile));
        spec.executeTest("testBreakBlocks", this);
    }

    @Test
    public void testSpanningDeletions() throws IOException {
        final String expectedFile = getToolTestDataDir() + "expected/" + "SpanningDeletions.vcf";
        IntegrationTestSpec  spec = new IntegrationTestSpec (
                " -o %s -R " + b37KGReference +
                        " -V " + getToolTestDataDir() + "spanningDel.1.g.vcf -V " + getToolTestDataDir() + "spanningDel.2.g.vcf",
                Arrays.asList(expectedFile));
        spec.executeTest("testSpanningDeletions", this);
    }

    @Test
    public void testMultipleSpanningDeletionsForOneSample() throws IOException {
        final String expectedFile = getToolTestDataDir() + "expected/" + "MultipleSpanningDeletionsForOneSample.vcf";
        IntegrationTestSpec  spec = new IntegrationTestSpec (
                " -o %s -R " + b37KGReference +
                        " -V " + getToolTestDataDir() + "spanningDel.many.g.vcf",
                Arrays.asList(expectedFile));
        spec. executeTest("testMultipleSpanningDeletionsForOneSample", this);
    }

    @Test
    public void testMultipleSpanningDeletionsForOneSampleHaploid() throws IOException {
        final String expectedFile = getToolTestDataDir() + "expected/" + "MultipleSpanningDeletionsForOneSampleHaploid.vcf";
        IntegrationTestSpec  spec = new IntegrationTestSpec (
                " -o %s -R " + b37KGReference +
                        " -V " + getToolTestDataDir() + "spanningDel.many.haploid.g.vcf",
                Arrays.asList(expectedFile));
        spec.executeTest("testMultipleSpanningDeletionsForOneSampleHaploid", this);
    }

    @Test
    public void testMultipleSpanningDeletionsForOneSampleTetraploid() throws IOException {
        final String expectedFile = getToolTestDataDir() + "expected/" + "MultipleSpanningDeletionsForOneSampleTetraploid.vcf";
        IntegrationTestSpec  spec = new IntegrationTestSpec (
                " -o %s -R " + b37KGReference +
                        " -V " + getToolTestDataDir() + "spanningDel.many.tetraploid.g.vcf",
                Arrays.asList(expectedFile));
                //Arrays.asList("0ec79471550ec5e30540f68cb0651b14"));
        spec.executeTest("testMultipleSpanningDeletionsForOneSampleTetraploid", this);
    }

    @Test
    public void testWrongReferenceBaseBugFix() throws Exception {
        final String expectedFile = getToolTestDataDir() + "expected/" + "WrongReferenceBaseBugFix.vcf";
        final String cmd = " -R " + b37KGReference + " -V " + (getToolTestDataDir() + "combine-gvcf-wrong-ref-input1.vcf"
                + " -V " + (getToolTestDataDir() + "combine-gvcf-wrong-ref-input2.vcf") + " -o %s ");
        final IntegrationTestSpec  spec = new IntegrationTestSpec (cmd, Arrays.asList(expectedFile));
        spec.executeTest("testWrongReferenceBaseBugFix",this);

    }

    @Test
    public void testBasepairResolutionInput() throws Exception {
        final String expectedFile = getToolTestDataDir() + "expected/" + "BasepairResolutionInput.vcf";
        final String cmd = " -R " + b37KGReference + " -o %s -V " + getToolTestDataDir() + "gvcf.basepairResolution.vcf";
        final IntegrationTestSpec spec = new IntegrationTestSpec (cmd, Arrays.asList(expectedFile));
        spec.executeTest("testBasepairResolutionInput", this);
    }
}
