package org.exbio.tfprio.steps.ehmm;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.File;
import java.util.*;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class ApplyModel extends ExecutableStep<Configs> {
    public final Map<String, OutputFile> outputDirs = new HashMap<>();
    private final Map<String, InputFile> regions = new HashMap<>();
    private final InputFile chromosomeSizes;
    private final InputFile constructedModel;
    private final RequiredConfig<File> bamDirectory = new RequiredConfig<>(configs.ehmm.bamDirectory);
    private final InputFile applyModelScript;

    public ApplyModel(Configs configs, Map<String, OutputFile> regions, OutputFile chromosomeSizes,
                      OutputFile constructedModel) {
        super(configs, false, regions.values(), chromosomeSizes, constructedModel);
        this.chromosomeSizes = addInput(chromosomeSizes);
        this.constructedModel = addInput(constructedModel);
        this.applyModelScript = addInput(getClass().getResourceAsStream("applyModel.R"), "applyModel.R");
        // apply model for each group and collapse bed files within
        regions.forEach((s, r) -> {
            this.regions.put(s, addInput(r));
            this.outputDirs.put(s, addOutput(s));
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                regions.forEach((s, r) -> add(() -> {
                    // build promoter model
                    String ehmmCommand = String.join(" ","Rscript",
                            applyModelScript.getAbsolutePath(),
                            "-r", r.getAbsolutePath(),
                            "-g", chromosomeSizes.getParent(),
                            "-m", constructedModel.getAbsolutePath(),
                            "-b", bamDirectory.get().getAbsolutePath(),
                            "-o", outputDirs.get(s).getAbsolutePath());
                    executeAndWait(ehmmCommand, false);
                    return true;
                }));
                return true;
            });
        }};
    }
}
