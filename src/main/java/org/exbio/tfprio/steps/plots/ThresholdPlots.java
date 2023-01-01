package org.exbio.tfprio.steps.plots;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class ThresholdPlots extends ExecutableStep {
    // Originally, there wer also plots for same stages and different stages here
    // Removed them, because it is not trivial to define what is a same stage and what is a different stage
    // First letter works just for the demo mouse data

    public final Map<String, Map<String, Map<Double, OutputFile>>> outputFiles = new HashMap<>();
    private final Map<InputFile, OutputFile> bridge = new HashMap<>();
    private final InputFile script;

    private final Map<InputFile, String> inputFileHm = new HashMap<>();
    private final Map<InputFile, String> inputFilePairing = new HashMap<>();

    public ThresholdPlots(Map<String, Map<String, Map<Double, OutputFile>>> coefficientFiles) {
        super(false, coefficientFiles.values().stream().flatMap(x -> x.values().stream()).flatMap(
                col -> col.values().stream()).toList());

        script = addInput(getClass().getResourceAsStream("thresholdPlot.py"), "thresholdPlot.py");

        coefficientFiles.forEach((pairing, hmMap) -> {
            OutputFile inputPairing = new OutputFile(inputDirectory, pairing);
            OutputFile outputPairing = new OutputFile(outputDirectory, pairing);

            hmMap.forEach((hm, thresholdMap) -> {
                OutputFile inputHm = new OutputFile(inputPairing, hm);
                OutputFile outputHm = new OutputFile(outputPairing, hm);

                thresholdMap.forEach((threshold, input) -> {
                    OutputFile outputFile = addOutput(outputHm, threshold + ".png");
                    InputFile inputFile = addInput(inputHm, input);
                    bridge.put(inputFile, outputFile);
                    outputFiles.computeIfAbsent(pairing, s -> new HashMap<>()).computeIfAbsent(hm,
                            s -> new HashMap<>()).put(threshold, outputFile);
                    inputFileHm.put(inputFile, hm);
                    inputFilePairing.put(inputFile, pairing);
                });
            });
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((input, output) -> add(() -> {
                String[] pairingSplit = inputFilePairing.get(input).split("_");
                String hm = inputFileHm.get(input);
                String group1 = pairingSplit[0];
                String group2 = pairingSplit[1];

                String command =
                        String.format("python3 %s %s %s %s %s %s", script.getAbsolutePath(), input.getAbsolutePath(),
                                output.getAbsolutePath(), hm, group1, group2);

                executeAndWait(command, true);
                return true;
            }));
        }};
    }
}
