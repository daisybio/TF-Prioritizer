package org.exbio.tfprio.steps.ehmm;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class CombinePredictions extends ExecutableStep<Configs> {
    public final Map<String, OutputFile> outputFiles = new HashMap<>();
    public final Map<String, OutputFile> stats = new HashMap<>();
    private final InputFile predictionDir;
    private final InputFile script;

    public CombinePredictions(Configs configs, Map<String, OutputFile> predictionParentDir) {
        super(configs, false, predictionParentDir.values());
        predictionDir = addInput(new OutputFile(predictionParentDir.values().iterator().next().getParent()));
        this.script = addInput(getClass().getResourceAsStream("combinePredictions.R"), "combinePredictions.R");
        this.outputFiles.put("enhancers", addOutput("predictedEnhancers.bed"));
        this.outputFiles.put("promoters", addOutput("predictedPromoters.bed"));
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
