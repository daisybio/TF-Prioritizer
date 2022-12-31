package org.exbio.tfprio.steps.Dynamite;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;

import java.util.*;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class PrepareForClassification extends ExecutableStep {
    public final Map<String, Map<String, OutputFile>> outputFiles = new HashMap<>();
    private final InputFile script;
    private final Map<InputFile, OutputFile> bridge = new HashMap<>();

    public PrepareForClassification(Map<String, Map<String, OutputFile>> integratedData) {
        super(false, integratedData.values().stream().flatMap(
                stringCollectionMap -> stringCollectionMap.values().stream()).toList());
        script = addInput(getClass().getResourceAsStream(
                        "/org/exbio/tfprio/steps/TEPIC/DYNAMITE/Scripts/prepareForClassification.R"),
                "prepareForClassification.R");
        integratedData.forEach((pairing, hmMap) -> {
            OutputFile inputPairing = new OutputFile(inputDirectory, pairing);
            OutputFile outputPairing = new OutputFile(outputDirectory, pairing);

            hmMap.forEach((hm, input) -> {
                InputFile inputFile = addInput(inputPairing, input);

                OutputFile outputDir = addOutput(outputPairing, hm);
                OutputFile outputFile = new OutputFile(outputDir, hm + ".tsv");
                bridge.put(inputFile, outputFile);
                outputFiles.computeIfAbsent(pairing, s -> new HashMap<>()).put(hm, outputDir);
            });
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((inputFile, outputFile) -> add(() -> {
                List<String> commandElements = List.of("Rscript", script.getAbsolutePath(), inputFile.getAbsolutePath(),
                        outputFile.getAbsolutePath());

                String command = String.join(" ", commandElements);

                executeAndWait(command, true);

                return true;
            }));
        }};
    }
}
