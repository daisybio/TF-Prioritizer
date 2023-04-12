package org.exbio.tfprio.steps.ehmm;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.util.Collection;
import java.util.HashSet;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class SplitLearningData extends ExecutableStep<Configs> {
    private final OutputFile outputDirectory = addOutput("out");
    public final OutputFile enhancers = addOutput(outputDirectory, "enhancers.bed");
    public final OutputFile promoters = addOutput(outputDirectory, "promoters.bed");
    public final OutputFile background = addOutput(outputDirectory, "background.bed");

    private final InputFile allEnhancers;
    private final InputFile allPromoters;
    private final InputFile allBackground;
    private final InputFile script;

    public SplitLearningData(Configs configs, OutputFile allEnhancers, OutputFile allPromoters,
                             OutputFile allBackground) {
        super(configs, false, allEnhancers, allPromoters, allBackground);
        this.allEnhancers = addInput(allEnhancers);
        this.allPromoters = addInput(allPromoters);
        this.allBackground = addInput(allBackground);
        this.script = addInput(getClass().getResourceAsStream("SplitTrainingData.R"), "SplitTrainingData.R");
    }


    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                // split learning data to background, enhancer, and promoter regions
                String ehmmCommand = String.join(" ","Rscript",
                        script.getAbsolutePath(),
                        "-bg", allBackground.getAbsolutePath(),
                        "-e", allEnhancers.getAbsolutePath(),
                        "-p", allPromoters.getAbsolutePath(),
                        "-o", outputDirectory.getAbsolutePath());
                executeAndWait(ehmmCommand, true);
                return true;
            });
        }};
    }
}
