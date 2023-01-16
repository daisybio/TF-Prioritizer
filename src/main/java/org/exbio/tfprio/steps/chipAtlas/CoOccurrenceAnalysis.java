package org.exbio.tfprio.steps.chipAtlas;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;

import java.util.Collection;
import java.util.HashSet;
import java.util.concurrent.Callable;

public class CoOccurrenceAnalysis extends ExecutableStep {
    private final InputFile inputDirectory;

    public CoOccurrenceAnalysis(OutputFile chipAtlasDirectory) {
        super(false, chipAtlasDirectory);

        inputDirectory = addInput(chipAtlasDirectory);
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                logger.warn("Not implemented yet");
                return true;
            });
        }};
    }
}
