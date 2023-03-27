package org.exbio.tfprio.steps.ehmm;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.util.Collection;
import java.util.HashSet;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;


public class ConstructModel extends ExecutableStep<Configs> {
    private final OutputFile outputDir = addOutput("constructedModel");
    public final OutputFile outputFile = addOutput(outputDir, "model.txt");
    private final InputFile backgroundModel;
    private final InputFile enhancerModel;
    private final InputFile promoterModel;

    private final RequiredConfig<Integer> nThreads = new RequiredConfig<>(configs.ehmm.nThreads);
    private final InputFile constructModelRscript;

    public ConstructModel(Configs configs, OutputFile backgroundModel, OutputFile enhancerModel, OutputFile promoterModel){
        super(configs, false, backgroundModel, enhancerModel, promoterModel);
        this.backgroundModel = addInput(backgroundModel);
        this.enhancerModel = addInput(enhancerModel);
        this.promoterModel = addInput(promoterModel);
        this.constructModelRscript = addInput(getClass().getResourceAsStream("constructModel.R"), "constructModel.R");
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                String ehmmCommand = String.join(" ","Rscript",
                        constructModelRscript.getAbsolutePath(),
                        "-b", backgroundModel.getAbsolutePath(),
                        "-e", enhancerModel.getAbsolutePath(),
                        "-p", promoterModel.getAbsolutePath(),
                        "-t", nThreads.toString(),
                        "-o", outputDir.getAbsolutePath());
                executeAndWait(ehmmCommand, false);
                return true;
            });
        }};
    }
}
