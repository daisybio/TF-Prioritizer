package org.exbio.tfprio.steps.distributionAnalysis;

import org.apache.commons.lang3.tuple.Pair;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class CreatePlots extends ExecutableStep<Configs> {
    // TODO: remove index from being stored
    public final Map<String, OutputFile> statFiles = new HashMap<>();
    public final Map<String, OutputFile> plotFiles = new HashMap<>();
    private final Map<Pair<InputFile, InputFile>, Pair<OutputFile, OutputFile>> bridge = new HashMap<>();

    private final InputFile script;

    public CreatePlots(Configs configs, Map<String, OutputFile> hmTfDirectory, Map<String, OutputFile> hmBackground) {
        super(configs, false, hmTfDirectory.values(), hmBackground.values());

        script = addInput(getClass().getResourceAsStream("plots.py"), "plots.py");

        OutputFile inTf = new OutputFile(inputDirectory, "tf");
        OutputFile inBackground = new OutputFile(inputDirectory, "background");

        OutputFile outStats = addOutput("stats");

        hmTfDirectory.forEach((hm, directory) -> {
            InputFile inputTf = addInput(inTf, directory);
            InputFile inputBackground = addInput(inBackground, hmBackground.get(hm));

            OutputFile outputHm = addOutput(hm);
            OutputFile outputStats = addOutput(outStats, hm + ".tsv");
            statFiles.put(hm, outputStats);
            plotFiles.put(hm, outputHm);

            bridge.put(Pair.of(inputTf, inputBackground), Pair.of(outputHm, outputStats));
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((input, output) -> add(() -> {
                InputFile inputTf = input.getLeft();
                InputFile inputBackground = input.getRight();
                OutputFile outputHm = output.getLeft();
                OutputFile outputStats = output.getRight();

                String command = "python3 " + script + " --backgroundFile " + inputBackground.getAbsolutePath() +
                        " --statsFile " + outputStats.getAbsolutePath() + " --inputDirectory " +
                        inputTf.getAbsolutePath() + " --outputDirectory " + outputHm.getAbsolutePath();

                executeAndWait(command, true);

                return true;
            }));
        }};
    }
}
