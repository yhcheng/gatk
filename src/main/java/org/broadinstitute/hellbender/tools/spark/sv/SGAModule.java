package org.broadinstitute.hellbender.tools.spark.sv;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Represents an SGA module that can be called via "run" (in base) to do actual work.
 */
public final class SGAModule extends CMDLineProgramModule{

    private static String programName = "sga";

    private final String moduleName;

    public SGAModule(final String moduleName){
        this.moduleName = moduleName;
    }

    @Override
    public ArrayList<String> initializeCommands() throws IOException {
        checkIfProgramIsAvailableOnHost(programName);
        final ArrayList<String> res = new ArrayList<>();
        res.add(programName);
        res.add(moduleName);
        return res;
    }

    void setupWorkingEnvironment(final ProcessBuilder builder, final String... args) {

    }
}
