package org.broadinstitute.hellbender.tools.spark.sv;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

/**
 * Tests RunMinimalBWAMEM.
 */
public final class RunMinimalBWAMEMTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "spark/sv");

    @Test(groups="sv")
    public void testSeparate() throws IOException {


    }

    @Test(groups="sv")
    public void testInterLeaved() throws IOException {

        try{
            CMDLineProgramModule.checkIfProgramIsAvailableOnHost("bwa");
        } catch(final IOException e){
            System.err.println(e.getMessage());
            return;
        }

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("-p");

        // IO arguments
        final File input = new File(TEST_DATA_DIR, "RunMinimalBWAMEMTest.fastq");
        args.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        args.add(input.getAbsolutePath());

        final File wkDir = BaseTest.createTempDir("dummy");
        args.add("-" + "outDir");
        args.add(wkDir.getAbsolutePath());

        final File output = new File(wkDir, "test.sam");
        output.createNewFile();
        args.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        args.add(output.getName());

        final File REF = new File(b37_reference_20_21);
        args.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
        args.add(REF.getAbsolutePath());

        this.runCommandLine(args.getArgsArray());
    }
}
