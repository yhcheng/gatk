package org.broadinstitute.hellbender.tools.spark.sv;


import java.util.ArrayList;
import java.io.IOException;

/**
 * Represents an BWA module that can be called via "run" (in base) to do actual alignment work.
 */
public final class BWAMEMModule extends CMDLineProgramModule{

    @Override
    public ArrayList<String> initializeCommands() throws IOException{
        checkIfProgramIsAvailableOnHost("bwa");
        final ArrayList<String> res = new ArrayList<>();
        res.add("bwa");
        res.add("mem");
        return res;
    }

    @Override
    public void setupWorkingEnvironment(final ProcessBuilder builder, final String... args){

    }
}