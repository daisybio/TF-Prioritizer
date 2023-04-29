package org.exbio.tfprio.steps.ehmm;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.OptionalConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.util.*;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class CombinePredictions extends ExecutableStep<Configs> {
    public final Map<String, OutputFile> outputFiles = new HashMap<>();
    public final Map<String, OutputFile> stats = new HashMap<>();

    private final OptionalConfig<Double> score = new OptionalConfig<>(configs.ehmm.score, false);
    private final InputFile predictionDir;
    private final InputFile script;

    public CombinePredictions(Configs configs, OutputFile predictionParentDir) {
        super(configs, false, predictionParentDir);
        predictionDir = addInput(predictionParentDir);
        this.script = addInput(getClass().getResourceAsStream("combinePredictions.R"), "combinePredictions.R");
        this.outputFiles.put("enhancers", addOutput("predictedEnhancers.bed"));
        this.outputFiles.put("promoters", addOutput("predictedPromoters.bed"));
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                List<String> cmd = new ArrayList<>(List.of("Rscript",
                        script.getAbsolutePath(),
                        "-d", predictionDir.getAbsolutePath(),
                        "-o", outputDirectory.getAbsolutePath()));
                if (score.isSet()) {
                    cmd.add("-t");
                    cmd.add(score.toString());
                }
                executeAndWait(cmd, false);
                return true;
            });
        }};
    }
}
