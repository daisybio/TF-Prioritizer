package org.exbio.tfprio.steps.ehmm;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.Callable;

public class CombinePredictions extends ExecutableStep<Configs> {
    public final OutputFile promoters = addOutput("predictedPromoters.bed");
    public final OutputFile enhancers = addOutput("predictedEnhancers.bed");
    private final Map<String, InputFile> predictionDirs = new HashMap<>();
    private final InputFile script;

    public CombinePredictions(Configs configs, Map<String, OutputFile> predictionDirs) {
        super(configs, false, predictionDirs.values());
        predictionDirs.forEach((cre, dir) -> {
            this.predictionDirs.put(cre, addInput(dir));
        });
        this.script = addInput(getClass().getResourceAsStream("combinePredictions.R"), "combinePredictions.R");
    }


    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return null;
    }
}
