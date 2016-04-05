package org.broadinstitute.hellbender.tools.spark.sv;


import java.util.ArrayList;
import java.util.List;
import java.io.IOException;

/**
 * Represents an BWA module that can be called via "run" (in base) to do actual alignment work.
 */
public final class BWAMEMModule extends CMDLineProgramModule{

    @Override
    public List<String> initializeCommands() throws IOException{
        checkIfProgramIsAvailableOnHost("bwa");
        final List<String> res = new ArrayList<>();
        res.add("bwa");
        res.add("mem");
        return res;
    }
}