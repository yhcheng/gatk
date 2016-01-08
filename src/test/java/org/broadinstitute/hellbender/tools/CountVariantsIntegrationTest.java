package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.CountVariants;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public final class CountVariantsIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return CountVariants.class.getSimpleName();
    }

    @Test(dataProvider = "countVariantsVCFInputs")
    public void testCountVariants(String fileIn, String moreArgs, long expectedCount) throws Exception {
        final File ORIG_FILE = new File(getTestDataDir(), fileIn);
        ArgumentsBuilder ab = new ArgumentsBuilder();
        ab.add("--variant " + ORIG_FILE.getAbsolutePath());
        ab.add(moreArgs);
        final Object res = this.runCommandLine(ab.getArgsArray());
        Assert.assertEquals(res, expectedCount);
    }

    @DataProvider(name="countVariantsVCFInputs")
    public Object[][] countVariantsVCFInputs() {
        return new Object[][]{
                {"count_variants.vcf", "", 26L},
                {"count_variants.vcf.blockgz.gz", "", 26L},
                {"count_variants_withSequenceDict.vcf", "", 26L},
                {"count_variants_withSequenceDict.vcf", "-L 1", 14L},
        };
    }

    @Test(expectedExceptions = UserException.class)
    public void testCountVariants_requiresIndexForInterval() throws Exception {
        String fileIn = "count_variants_withSequenceDict_noIndex.vcf";
        String moreArgs = "-L 1";
        final File ORIG_FILE = new File(getTestDataDir(), fileIn);
        ArgumentsBuilder ab = new ArgumentsBuilder();
        ab.add("--variant " + ORIG_FILE.getAbsolutePath());
        ab.add(moreArgs);
        final Object res = this.runCommandLine(ab.getArgsArray());
    }

    @Test(expectedExceptions = UserException.class)
    public void testCountVariants_missingContigForInterval() throws Exception {
        String fileIn = "count_variants_withSequenceDict.vcf";
        String moreArgs = "-L 25";
        final File ORIG_FILE = new File(getTestDataDir(), fileIn);
        ArgumentsBuilder ab = new ArgumentsBuilder();
        ab.add("--variant " + ORIG_FILE.getAbsolutePath());
        ab.add(moreArgs);
        final Object res = this.runCommandLine(ab.getArgsArray());
    }

    @Test(expectedExceptions = UserException.class)
    public void testCountVariants_requiresSequenceDictionaryForInterval() throws Exception {
        String fileIn = "count_variants.vcf";
        String moreArgs = "-L 1";
        final File ORIG_FILE = new File(getTestDataDir(), fileIn);
        ArgumentsBuilder ab = new ArgumentsBuilder();
        ab.add("--variant " + ORIG_FILE.getAbsolutePath());
        ab.add(moreArgs);
        final Object res = this.runCommandLine(ab.getArgsArray());
    }

}
