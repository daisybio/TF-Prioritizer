package org.exbio.tfprio.steps.ehmm;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.util.Collection;
import java.util.HashSet;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class CombinePredictions extends ExecutableStep<Configs> {
    public final OutputFile promoters = addOutput("predictedPromoters.bed");
    public final OutputFile enhancers = addOutput("predictedEnhancers.bed");
    private final InputFile predictionDir;
    private final InputFile script;

    public CombinePredictions(Configs configs, OutputFile predictionParentDir) {
        super(configs, false, predictionParentDir);
        predictionDir = addInput(predictionParentDir);
        this.script = addInput(getClass().getResourceAsStream("combinePredictions.R"), "combinePredictions.R");
    }


    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                String cmd = String.join(" ","Rscript",
                        script.getAbsolutePath(),
                        "-d", predictionDir.getAbsolutePath(),
                        "-o", outputDirectory.getAbsolutePath());
                executeAndWait(cmd, false);
                return true;
            });
        }};
    }
}
