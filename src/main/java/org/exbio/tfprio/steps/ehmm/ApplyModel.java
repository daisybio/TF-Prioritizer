package org.exbio.tfprio.steps.ehmm;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.File;
import java.util.*;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.FileManagement.extend;
import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class ApplyModel extends ExecutableStep<Configs> {
    public final Map<String, OutputFile> outputDirs = new HashMap<>();
    private final Map<String, InputFile> regions = new HashMap<>();
    private final Map<String, InputFile> bamDirs = new HashMap<>();
    private final InputFile chromosomeSizes;
    private final InputFile constructedModel;
    private final RequiredConfig<File> bamDirectory = new RequiredConfig<>(configs.ehmm.bamDirectory);
    private final RequiredConfig<Integer> nThreads = new RequiredConfig<>(configs.ehmm.nThreads);
    private final RequiredConfig<Double> pseudoCounts = new RequiredConfig<>(configs.ehmm.pseudoCount);
    private final Integer binSize = new RequiredConfig<>(configs.ehmm.nBins).get();
    private final InputFile applyModelScript;

    public ApplyModel(Configs configs, Map<String, OutputFile> regions, OutputFile chromosomeSizes,
                      OutputFile constructedModel) {
        super(configs, false, regions.values(), chromosomeSizes, constructedModel);
        this.chromosomeSizes = addInput(chromosomeSizes);
        this.constructedModel = addInput(constructedModel);
        this.applyModelScript = addInput(getClass().getResourceAsStream("applyModel.R"), "applyModel.R");

        regions.forEach((s, r) -> {
            OutputFile groupDir = new OutputFile(inputDirectory, s);
            this.regions.put(s, addInput(groupDir, r));
            this.bamDirs.put(s, addInput(groupDir, new OutputFile(
                    extend(bamDirectory.get(), s).getAbsolutePath())));
            this.outputDirs.put(s, addOutput(s));
        });
        Arrays.stream(Objects.requireNonNull(this.bamDirectory.get().listFiles()))
                .filter(File::isDirectory)
                .forEach(sDir -> {
                    OutputFile groupIn = new OutputFile(inputDirectory, sDir.getName());
                    Arrays.stream(Objects.requireNonNull(sDir.listFiles()))
                            .filter(File::isDirectory)
                            .forEach(hmDir -> {
                                    OutputFile hmIn = new OutputFile(groupIn, hmDir.getName());
                                    Arrays.stream(Objects.requireNonNull(hmDir.listFiles()))
                                            .filter(f -> f.toPath().endsWith(".bam"))
                                            .forEach(f -> {
                                                OutputFile bamFile = new OutputFile(f.getAbsolutePath());
                                                OutputFile baiFile = new OutputFile(bamFile.getAbsolutePath()+".bai");
                                                addInput(hmIn, bamFile);
                                                if (baiFile.exists()) addInput(hmIn, baiFile);
                                            });
                            });
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
