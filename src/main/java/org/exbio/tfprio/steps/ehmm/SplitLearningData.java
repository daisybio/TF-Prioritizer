package org.exbio.tfprio.steps.ehmm;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.File;
import java.util.Collection;
import java.util.HashSet;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class SplitLearningData extends ExecutableStep<Configs> {
    private final OutputFile outputDirectory = addOutput("out");
    public final OutputFile trainEnhancers = addOutput(outputDirectory, "trainEnhancers.bed");
    public final OutputFile trainPromoters = addOutput(outputDirectory, "trainPromoters.bed");
    public final OutputFile trainBackground = addOutput(outputDirectory, "trainBackground.bed");
    public final OutputFile testEnhancers = addOutput(outputDirectory, "testEnhancers.bed");
    public final OutputFile testPromoters = addOutput(outputDirectory, "testPromoters.bed");
    public final OutputFile testBackground = addOutput(outputDirectory, "testBackground.bed");

    private final InputFile allEnhancers;
    private final InputFile allPromoters;
    private final InputFile allBackground;
    private final RequiredConfig<File> gtfFile = new RequiredConfig<>(configs.inputConfigs.geneAnnotationFile);
    private final RequiredConfig<Integer> genomicRegionSize = new RequiredConfig<>(configs.ehmm.genomicRegionSize);
    private final RequiredConfig<Double> trainSplit = new RequiredConfig<>(configs.ehmm.trainSplit);
    private final InputFile gtf;
    private final InputFile script;

    public SplitLearningData(Configs configs, OutputFile allEnhancers, OutputFile allPromoters,
                             OutputFile allBackground) {
        super(configs, false, allEnhancers, allPromoters, allBackground);
        this.allEnhancers = addInput(allEnhancers);
        this.allPromoters = addInput(allPromoters);
        this.allBackground = addInput(allBackground);
        this.gtf = addInput(gtfFile);
        this.script = addInput(getClass().getResourceAsStream("SplitTrainingData.R"), "SplitTrainingData.R");
    }


    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                // split learning data to background, enhancer, and promoter regions
                String ehmmCommand = String.join(" ","Rscript",
                        script.getAbsolutePath(),
                        "-b", allBackground.getAbsolutePath(),
                        "-e", allEnhancers.getAbsolutePath(),
                        "-p", allPromoters.getAbsolutePath(),
                        "-g", gtf.getAbsolutePath(),
                        "-n", trainSplit.toString(),
                        "-s", genomicRegionSize.toString(),
                        "-o", outputDirectory.getAbsolutePath());
                executeAndWait(ehmmCommand, true);
                return true;
            });
        }};
    }
}
