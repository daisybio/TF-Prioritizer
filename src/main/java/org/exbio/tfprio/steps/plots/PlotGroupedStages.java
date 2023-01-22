package org.exbio.tfprio.steps.plots;

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
import java.util.stream.Stream;

import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class PlotGroupedStages extends ExecutableStep<Configs> {
    public final Map<String, Map<Double, Pair<OutputFile, OutputFile>>> outputFiles = new HashMap<>();
    private final Map<InputFile, OutputFile> bridge = new HashMap<>();
    private final InputFile script;

    public PlotGroupedStages(Configs configs, Map<String, Map<Double, Pair<OutputFile, OutputFile>>> groupedStages) {
        super(configs, false, groupedStages.values().stream().flatMap(x -> x.values().stream()).flatMap(
                pair -> Stream.of(pair.getLeft(), pair.getRight())).toList());

        script = addInput(getClass().getResourceAsStream("groupedStagePlot.py"), "groupedStagePlot.py");

        groupedStages.forEach((hm, thresholdMap) -> {
            OutputFile inHm = new OutputFile(inputDirectory, hm);
            OutputFile outHm = new OutputFile(outputDirectory, hm);
            thresholdMap.forEach((threshold, pair) -> {
                OutputFile inThreshold = new OutputFile(inHm, threshold.toString());
                OutputFile outThreshold = new OutputFile(outHm, threshold.toString());

                InputFile inputSame = addInput(inThreshold, pair.getLeft());
                InputFile inputDifferent = addInput(inThreshold, pair.getRight());

                OutputFile outputSame = addOutput(outThreshold, "same.png");
                OutputFile outputDifferent = addOutput(outThreshold, "different.png");

                outputFiles.computeIfAbsent(hm, s -> new HashMap<>()).put(threshold,
                        Pair.of(outputSame, outputDifferent));

                bridge.put(inputSame, outputSame);
                bridge.put(inputDifferent, outputDifferent);
            });
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((input, output) -> add(() -> {
                String command = "python3 " + script + " " + input + " " + output;
                executeAndWait(command, true);
                return true;
            }));
        }};
    }
}
