package org.exbio.tfprio.steps.ehmm;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.File;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class ApplyModel extends ExecutableStep<Configs> {
    public final Map<String, OutputFile> outputDirs = new HashMap<>();
    private final Map<String, InputFile> regions = new HashMap<>();
    private final InputFile chromosomeSizes;
    private final InputFile constructedModel;
    private final Map<String, InputFile> bamDirs;
    private final RequiredConfig<Integer> nThreads = new RequiredConfig<>(configs.ehmm.nThreads);
    private final RequiredConfig<Double> pseudoCounts = new RequiredConfig<>(configs.ehmm.pseudoCount);
    private final Integer binSize = new RequiredConfig<>(configs.ehmm.nBins).get();
    private final InputFile applyModelScript;

    public ApplyModel(Configs configs, Map<String, OutputFile> regions, OutputFile chromosomeSizes,
                      OutputFile constructedModel, Map<String, OutputFile> bamDirs) {
        super(configs, false, Stream.of(regions, bamDirs)
                .flatMap(m -> m.values().stream()).collect(Collectors.toSet()), chromosomeSizes, constructedModel);
        this.chromosomeSizes = addInput(chromosomeSizes);
        this.constructedModel = addInput(constructedModel);
        this.bamDirs = bamDirs.entrySet().stream()
                .map(e -> new AbstractMap
                        .SimpleEntry<>(e.getKey(), addInput(e.getValue())))
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
        this.applyModelScript = addInput(getClass().getResourceAsStream("applyModel.R"), "applyModel.R");
        regions.forEach((s, r) -> {
            OutputFile groupDir = new OutputFile(inputDirectory, s);
            this.regions.put(s, addInput(groupDir, r));
            this.outputDirs.put(s, addOutput(s));
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            regions.forEach((s, r) -> add(() -> {
                // apply model to each group
                String ehmmCommand = String.join(" ","Rscript",
                        applyModelScript.getAbsolutePath(),
                        "-r", r.getAbsolutePath(),
                        "-g", chromosomeSizes.getAbsolutePath(),
                        "-m", constructedModel.getAbsolutePath(),
                        "-b", bamDirs.get(s).getAbsolutePath(),
                        "-o", outputDirs.get(s).getAbsolutePath(),
                        "-s", binSize.toString(),
                        "-t", nThreads.toString(),
                        "-p", pseudoCounts.toString());
                executeAndWait(ehmmCommand, false);
                return true;
            }));
        }};
    }
}
