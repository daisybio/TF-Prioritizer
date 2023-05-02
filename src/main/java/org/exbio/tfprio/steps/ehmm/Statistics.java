package org.exbio.tfprio.steps.ehmm;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.util.*;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class Statistics extends ExecutableStep<Configs> {
    public final OutputFile outputFile = addOutput("Statistics.tsv");

    private final InputFile collapsedPeakParents;

    private final InputFile enhancerPrediction;
    private final InputFile enhancerReference;

    private final InputFile promoterPrediction;
    private final InputFile promoterReference;

    private final InputFile script;

    public Statistics(Configs configs, Map<String, OutputFile> combinedPredictions, OutputFile collapsedPeakParents,
                      OutputFile enhancerReference, OutputFile promoterReference) {
        super(configs, false, combinedPredictions.values(), collapsedPeakParents);

        this.collapsedPeakParents = addInput(collapsedPeakParents);
        this.script = addInput(getClass().getResourceAsStream("stats.R"), "stats.R");
        this.enhancerPrediction = addInput(combinedPredictions.get("enhancers"));
        this.promoterPrediction = addInput(combinedPredictions.get("promoters"));
        this.enhancerReference = addInput(enhancerReference);
        this.promoterReference = addInput(promoterReference);
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                List<String> cmd = new ArrayList<>(List.of("Rscript",
                        script.getAbsolutePath(),
                        "--p_enhancers", enhancerPrediction.getAbsolutePath(),
                        "--r_enhancers", enhancerReference.getAbsolutePath(),
                        "--p_promoters", promoterPrediction.getAbsolutePath(),
                        "--r_promoters", promoterReference.getAbsolutePath(),
                        "-i", collapsedPeakParents.getAbsolutePath(),
                        "-o", outputFile.getAbsolutePath()));
                executeAndWait(cmd, false);
                return true;
            });
        }};
    }
}
