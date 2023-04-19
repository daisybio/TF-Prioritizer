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
    private final InputFile chromosomeSizes;
    private final InputFile constructedModel;
    private final InputFile bamDirectory;
    private final RequiredConfig<Integer> nThreads = new RequiredConfig<>(configs.ehmm.nThreads);
    private final RequiredConfig<Double> pseudoCounts = new RequiredConfig<>(configs.ehmm.pseudoCount);
    private final Integer binSize = new RequiredConfig<>(configs.ehmm.nBins).get();
    private final InputFile applyModelScript;

    public ApplyModel(Configs configs, Map<String, OutputFile> regions, OutputFile chromosomeSizes,
                      OutputFile constructedModel) {
        super(configs, false, regions.values(), chromosomeSizes, constructedModel);
        this.chromosomeSizes = addInput(chromosomeSizes);
        this.constructedModel = addInput(constructedModel);
        this.bamDirectory = addInput(new OutputFile(new RequiredConfig<>(configs.ehmm.bamDirectory).get().getAbsolutePath()));
        this.applyModelScript = addInput(getClass().getResourceAsStream("applyModel.R"), "applyModel.R");

        regions.forEach((s, r) -> {
            OutputFile groupDir = new OutputFile(inputDirectory, s);
            this.regions.put(s, addInput(groupDir, r));
            this.outputDirs.put(s, addOutput(s));
        });
        // add bam and bai files as inputs
        Arrays.stream(Objects.requireNonNull(bamDirectory.listFiles()))
                .forEach(group -> Arrays.stream(Objects.requireNonNull(bamDirectory.listFiles()))
                        .filter(f -> f.toString().endsWith(".bam"))
                        .forEach(b -> {
                            addInput(new OutputFile(b.getAbsolutePath()));
                            addInput(new OutputFile(b.getAbsolutePath()+".bai"));
                        }));
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            regions.forEach((s, r) -> add(() -> {
                File bamDir = extend(bamDirectory, s);
                // apply model to each group
                String ehmmCommand = String.join(" ","Rscript",
                        applyModelScript.getAbsolutePath(),
                        "-r", r.getAbsolutePath(),
                        "-g", chromosomeSizes.getAbsolutePath(),
                        "-m", constructedModel.getAbsolutePath(),
                        "-b", bamDir.getAbsolutePath(),
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
