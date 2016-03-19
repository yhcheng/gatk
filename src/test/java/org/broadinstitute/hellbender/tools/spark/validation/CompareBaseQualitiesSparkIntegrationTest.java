package org.broadinstitute.hellbender.tools.spark.validation;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public final class CompareBaseQualitiesSparkIntegrationTest extends CommandLineProgramTest {
    final String resourceDir = getTestDataDir() + "/validation/";
    final String singleReadDiffFileName = "single.read.qual.diff.txt";
    final String twoReadDiffFileName = "two.reads.qual.diff.txt";

    @Test
    public void identicalBamTest() throws Exception {
        // This test verifies that if two bams are identical that they produce zero diffs.
        final File firstBam = new File(resourceDir, "CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam");
        final File secondBam = new File(resourceDir, "CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam");

        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(firstBam.getCanonicalPath());
        args.add("--" + "I2");
        args.add(secondBam.getCanonicalPath());
        args.add("--" + "throwOnDiff");
        args.add("true");
        this.runCommandLine(args.getArgsArray());
    }

    @DataProvider(name = "CompareBasesProvider")
    public Object[][] makeCompareBasesProvider() {
        final File outFile = BaseTest.createTempFile(getTestDataDir() + "temp.diff", "txt");
        final File firstBam = new File(resourceDir, "single.read.bam");
        final File secondBam = new File(resourceDir, "another.single.read.bam");
        final File firstCram = new File(resourceDir, "single.read.cram");
        final File secondCram = new File(resourceDir, "another.single.read.cram");
        final File referenceFile = new File(b37_reference_20_21);
        return new Object[][]{
                {firstBam, firstBam, null, outFile, singleReadDiffFileName},
                {firstBam, secondBam, null, outFile, twoReadDiffFileName},
                {firstCram, secondCram, referenceFile, outFile, twoReadDiffFileName},
                {firstBam, secondCram, referenceFile, outFile, twoReadDiffFileName},
                {firstCram, secondBam, referenceFile, outFile, twoReadDiffFileName},
        };

    }

    private void doTest(File firstBam, File secondBam, File referenceFile, File outFile, String diffFile) throws IOException {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(firstBam.getCanonicalPath());
        args.add("--" + "I2");
        args.add(secondBam.getCanonicalPath());
        args.add("--" + "O");
        args.add(outFile.getCanonicalPath());
        if (null != referenceFile) {
            args.add("-R");
            args.add(referenceFile.getAbsolutePath());
        }

        this.runCommandLine(args.getArgsArray());
    }

    @Test(expectedExceptions = UserException.MissingReference.class)
    public void nonExistentReferenceTest() throws Exception {
        final File firstCram = new File(resourceDir, "single.read.cram");
        final File secondCram = new File(resourceDir, "another.single.read.cram");

        doTest(firstCram, secondCram,
                new File("Nonexistent.fasta"),
                BaseTest.createTempFile(getTestDataDir() + "temp.diff", "txt"),
                twoReadDiffFileName);
    }

    @Test(dataProvider = "CompareBasesProvider")
    public void singleReadDiffTest(File firstBam, File secondBam, File referenceFile, File outFile, String diffFile) throws Exception {
        doTest(firstBam, secondBam, referenceFile,  outFile, diffFile);

        final File expected = new File(resourceDir, diffFile);
        IntegrationTestSpec.assertEqualTextFiles(outFile, expected);
    }

}