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

public class LearnPromoterModel extends ExecutableStep<Configs> {
    private final OutputFile outputDir = addOutput("out");
    public final OutputFile outputFile = addOutput(outputDir,"PromoterModel.RData");
    private final InputFile bedFile;
    private final InputFile bamDir;

    private final RequiredConfig<Integer> nStates = new RequiredConfig<>(configs.ehmm.nStates);
    private final RequiredConfig<Double> pseudoCount = new RequiredConfig<>(configs.ehmm.pseudoCount);
    private final RequiredConfig<Integer> nBins = new RequiredConfig<>(configs.ehmm.nBins);
    private final InputFile learnModelRscript;

    public LearnPromoterModel(Configs configs, OutputFile bedFile, OutputFile bamDir){
        super(configs, false, bedFile, bamDir);
        this.bedFile = addInput(bedFile);
        this.bamDir = addInput(bamDir);
        this.learnModelRscript = addInput(getClass().getResourceAsStream("learnModel.R"), "learnModel.R");
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                // build promoter model
                String ehmmCommand = String.join(" ","Rscript",
                        learnModelRscript.getAbsolutePath(),
                        "-r", bedFile.getAbsolutePath(),
                        "-m", bamDir.getAbsolutePath(),
                        "-f", "PromoterModel",
                        "-n", nStates.toString(),
                        "-b", nBins.toString(),
                        "-p", pseudoCount.toString(),
                        "-o", outputDir.getAbsolutePath());
                executeAndWait(ehmmCommand, true);
                return true;
            });
        }};
    }
}
