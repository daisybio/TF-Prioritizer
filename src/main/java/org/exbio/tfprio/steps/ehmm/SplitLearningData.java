package org.exbio.tfprio.steps.ehmm;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.OptionalConfig;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.File;
import java.util.*;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class SplitLearningData extends ExecutableStep<Configs> {
    private final OutputFile outputDirectory = addOutput("out");
    public final OutputFile trainEnhancers = addOutput(outputDirectory, "trainEnhancers.bed");
    public final OutputFile trainPromoters = addOutput(outputDirectory, "trainPromoters.bed");
    public final OutputFile trainBackground = addOutput(outputDirectory, "trainBackground.bed");

    private final InputFile allEnhancers;
    private final InputFile allPromoters;
    private final InputFile allBackground;
    private final RequiredConfig<File> gtfFile = new RequiredConfig<>(configs.inputConfigs.geneAnnotationFile);
    private final RequiredConfig<Integer> genomicRegionSize = new RequiredConfig<>(configs.ehmm.genomicRegionSize);
    private final OptionalConfig<Double> trainSplit = new OptionalConfig<>(configs.ehmm.trainSplit, false);
    private final OptionalConfig<Integer> nSamples = new OptionalConfig<>(configs.ehmm.nSamples, false);
    private final OptionalConfig<Integer> topQuantile = new OptionalConfig<>(configs.ehmm.topQuantile, false);
    private final OptionalConfig<Boolean> chipBG = new OptionalConfig<>(configs.ehmm.chipBG, false);
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
                List<String> ehmmCommand = new ArrayList<>(List.of("Rscript",
                        script.getAbsolutePath(),
                        "-b", allBackground.getAbsolutePath(),
                        "-e", allEnhancers.getAbsolutePath(),
                        "-p", allPromoters.getAbsolutePath(),
                        "-g", gtf.getAbsolutePath(),
                        "-s", genomicRegionSize.toString(),
                        "-o", outputDirectory.getAbsolutePath()));
                if (trainSplit.isSet()) {
                    ehmmCommand.add("-t");
                    ehmmCommand.add(trainSplit.toString());
                }
                if (nSamples.isSet()) {
                    ehmmCommand.add("-f");
                    ehmmCommand.add(nSamples.toString());
                }
                if (topQuantile.isSet()) {
                    ehmmCommand.add("-q");
                    ehmmCommand.add(topQuantile.toString());
                }
                if (chipBG.get()) {
                    ehmmCommand.add("--chip_bg");
                }
                executeAndWait(ehmmCommand, true);
                return true;
            });
        }};
    }
}
