package org.broadinstitute.hellbender.tools.spark.sv;

import org.apache.logging.log4j.LogManager;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;

import java.io.*;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;

// TODO: choose which parameters allowed to be tunable
// TODO: if throws, would temp files be cleaned up automatically?
// TODO: how to log processes' progression in Spark? Currently only information when things go seriously wrong is logged in exception messages.
@CommandLineProgramProperties(
        summary        = "Minimal program to call BWAMEM for performing alignment, allowing limited options.",
        oneLineSummary = "Minimal work to call BWAMEM.",
        programGroup   = StructuralVariationSparkProgramGroup.class)
public final class RunMinimalBWAMEM extends CommandLineProgram {

    @Argument(doc       = "Absolute path to input FASTQ/A file to be aligned.",
              shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
              fullName  = StandardArgumentDefinitions.INPUT_LONG_NAME,
              optional  = false)
    public String input = null;

    @Argument(doc       = "Absolute path to second input FASTQ/A file to be aligned." +
                            "Sequences in this file are assumed to have the same read name as they appear in the first input, and in the same order.",
              shortName = "I2",
              fullName  = "secondInput",
              optional  = true)
    public String secondInput = null;

    @Argument(doc       = "If set to true, performs smart pairing (input FASTQ assumed interleaved), and ignores second fastq input.",
              shortName = "p",
              fullName  = "interLeaved",
              optional  = true)
    public boolean interLeaved = false;

    @Argument(doc       = "Absolute path to reference of the target organism, if alignment of assembled contigs to reference is desired.",
              shortName = StandardArgumentDefinitions.REFERENCE_SHORT_NAME,
              fullName  = StandardArgumentDefinitions.REFERENCE_LONG_NAME,
              optional  = false)
    public String reference = null;

    @Argument(doc       = "Sam file to write results to.",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName  = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional  = false)
    public String samOutput = null;

    @Argument(doc       = "A path to a directory to write results to.",
              shortName = "outDir",
              fullName  = "outputDirectory",
              optional  = false)
    public String outputDir = null;

    @Argument(doc       = "Number of threads to use when running bwa.",
              shortName = "t",
              fullName  = "threads",
              optional  = true)
    public int threads = 1;

    @Argument(doc       = "Number of threads to use when running bwa.",
              shortName = "K",
              fullName  = "chunkSizeEachThread",
              optional  = true)
    public long chunkSize = 0;

    @Override
    public String doWork() throws RuntimeException{

        if(!interLeaved && secondInput==null){
            throw new RuntimeException("Missing cmd line argument options: \n" +
                                       "Option \"-p\" not set but only one FASTQ/A given.");
        }

        String stderrMessage = new String( LogManager.getLogger(this.getClass()).toString() );

        try{
            final BWAMEMModule bwamem = new BWAMEMModule();

            final ArrayList<String> bwaArgs = makeArgs();

            final File wkDir = new File(outputDir);
            final File samFile = new File(wkDir, samOutput);

            bwamem.run(bwaArgs, wkDir, samFile, stderrMessage);
        } catch(final IOException e){
            throw new RuntimeException(e.getMessage());
        } catch(final InterruptedException e){
            throw new RuntimeException(e.getMessage());
        }

        System.out.println(stderrMessage);
        return stderrMessage;
    }

    private ArrayList<String> makeArgs(){
        final ArrayList<String> args = new ArrayList<>();

        final Path pathToReference = Paths.get(reference);
        final Path pathToInput = Paths.get(input);

        args.add("t");
        args.add(Integer.toString(threads));

        if(0!=chunkSize){
            args.add("K");
            args.add(Long.toString(this.chunkSize));
        }

        if(interLeaved){ // paired reads, interleaved
            args.add("p");
            args.add(pathToReference.toString());
            args.add(pathToInput.toString());
        }else{
            if(null!=secondInput){ // paired reads, separate files
                final Path pathToSecondInput = Paths.get(secondInput);
                args.add(pathToSecondInput.toString());
            }else{       // SE reads
                args.add("S"); // skips mate rescuing and pairing
                args.add("P");
                args.add(pathToReference.toString());
                args.add(pathToInput.toString());
            }
        }
        return args;
    }
}