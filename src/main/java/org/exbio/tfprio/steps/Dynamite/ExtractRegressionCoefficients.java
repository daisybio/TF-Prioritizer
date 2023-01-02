package org.exbio.tfprio.steps.Dynamite;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;

import java.io.File;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.FileManagement.softLink;

public class ExtractRegressionCoefficients extends ExecutableStep {

    public final Map<String, Map<String, OutputFile>> outputFiles = new HashMap<>();
    private final Map<InputFile, OutputFile> bridge = new HashMap<>();

    public ExtractRegressionCoefficients(Map<String, Map<String, OutputFile>> rawDynamiteOutput) {
        super(false, rawDynamiteOutput.values().stream().flatMap(
                stringCollectionMap -> stringCollectionMap.values().stream()).toList());

        rawDynamiteOutput.forEach((pairing, hmMap) -> {
            OutputFile inputPairing = new OutputFile(inputDirectory, pairing);
            OutputFile outputPairing = new OutputFile(outputDirectory, pairing);

            hmMap.forEach((hm, inDir) -> {
                InputFile inputDir = addInput(inputPairing, inDir);
                OutputFile outputDir = addOutput(outputPairing, hm + ".tsv");

                bridge.put(inputDir, outputDir);
                outputFiles.computeIfAbsent(pairing, s -> new HashMap<>()).put(hm, outputDir);
            });
        });
    }

    @Override
    protected boolean doCreateFiles() {
        return false;
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((input, output) -> add(() -> {
                File[] considered = input.listFiles(file -> file.getName().equals(
                        "Regression_Coefficients_Entire_Data_Set_" + input.getName() + ".txt"));

                if (considered == null || considered.length == 0) {
                    return false;
                }

                File regressionCoefficients = considered[0];

                softLink(output, regressionCoefficients);

                return true;
            }));
        }};
    }
}
