package org.exbio.tfprio.steps.distributionAnalysis;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;

import java.util.Collection;
import java.util.HashSet;
import java.util.concurrent.Callable;

public class CreateHeatmaps extends ExecutableStep {
    private final InputFile script;

    // TODO: Implement this step
    // Store the normalized gene counts in the DESEQ step and re-use them here
    public CreateHeatmaps() {
        super();

        script = addInput(getClass().getResourceAsStream("heatmaps.R"), "heatmaps.R");
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
