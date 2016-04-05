package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.annotations.VisibleForTesting;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

public abstract class CMDLineProgramModule {

    public CMDLineProgramModule(){  }

    abstract List<String> initializeCommands() throws IOException;

    public void run(final String[] runtimeArguments, final File directoryToWorkIn, final File stdoutDestination) throws IOException, InterruptedException, RuntimeException{

        final List<String> commands = initializeCommands();
        for(final String arg : runtimeArguments){ commands.add(arg); }

        final ProcessBuilder builder = new ProcessBuilder(commands);
        setupEnvironment(builder, directoryToWorkIn, stdoutDestination);

        Process runProcess = builder.start();
        int exitStatus = runProcess.waitFor();

        if(0!=exitStatus){
            onError(commands, redirectSTDERR(runProcess), exitStatus);
        }
    }

    protected void setupEnvironment(final ProcessBuilder builder, final File directoryToWorkIn, final File stdoutDestination){
        // setup working directory
        // TODO: see if working environment needs to be set up
        builder.directory(directoryToWorkIn);
        if(null!=stdoutDestination) {
            builder.redirectOutput(stdoutDestination);
        }
    }

    protected static String redirectSTDERR(final Process process) throws IOException{
        final BufferedReader reader = new BufferedReader(new InputStreamReader(process.getErrorStream()));
        StringBuilder out = new StringBuilder();
        String line       = null;
        String previous   = null;
        while ((line = reader.readLine()) != null) {
            if (!line.equals(previous)) {
                previous = line;
                out.append(line).append('\n');
            }
        }
        return out.toString();
    }

    private static void onError(final List<String> commands, final String commandMessage, final int commandExitStatus) throws RuntimeException{
        String errorMessage = "";
        for(final String mess : commands){ errorMessage += " " + mess; }
        errorMessage += "\n" + commandMessage;
        throw new RuntimeException("Errors occurred while running: " + errorMessage +
                                    "\nWith exit status: " + commandExitStatus);
    }

    @VisibleForTesting
    static void checkIfProgramIsAvailableOnHost(final String programName) throws IOException{
        try{
            final List<String> commands = new ArrayList<>();
            commands.add("which");
            commands.add(programName);
            final ProcessBuilder testIfProgramIsAvailableOnHost = new ProcessBuilder(commands);
            final Process programPath = testIfProgramIsAvailableOnHost.start();
            int exitStatus = programPath.waitFor();
            if(0!=exitStatus){
                throw new IOException("Can't find " + programName + " programs in $PATH of host.");
            }
        } catch(final InterruptedException e){
            throw new IOException("Errors occurred while testing where " + programName + " lives. " + e.getMessage());
        }
    }
}
