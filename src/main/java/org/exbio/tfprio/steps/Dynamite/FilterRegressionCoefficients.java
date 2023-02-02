package org.exbio.tfprio.steps.Dynamite;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.*;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.Callable;

public class FilterRegressionCoefficients extends ExecutableStep<Configs> {
    public final Map<String, Map<String, OutputFile>> outputFiles = new HashMap<>();
    private final Map<InputFile, OutputFile> bridge = new HashMap<>();

    private final RequiredConfig<Double> minRegressionCoefficient =
            new RequiredConfig<>(configs.dynamite.minRegressionCoefficient);

    public FilterRegressionCoefficients(Configs configs, Map<String, Map<String, OutputFile>> pairingHmCoefficients) {
        super(configs, false, pairingHmCoefficients.values().stream().flatMap(
                stringCollectionMap -> stringCollectionMap.values().stream()).toList());

        pairingHmCoefficients.forEach((pairing, hmMap) -> {
            OutputFile inPairing = new OutputFile(inputDirectory, pairing);
            OutputFile outPairing = new OutputFile(outputDirectory, pairing);

            hmMap.forEach((hm, inFile) -> {
                InputFile input = addInput(inPairing, inFile);
                OutputFile outputDir = addOutput(outPairing, hm + ".tsv");

                bridge.put(input, outputDir);
                outputFiles.computeIfAbsent(pairing, s -> new HashMap<>()).put(hm, outputDir);
            });
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((input, output) -> add(() -> {
                try (var reader = new BufferedReader(new FileReader(input));
                     var writer = new BufferedWriter(new FileWriter(output))) {
                    writer.write(reader.readLine());

                    reader.lines().filter(line -> {
                        String[] split = line.split("\t");
                        double regressionCoefficient = Double.parseDouble(split[1]);
                        return Math.abs(regressionCoefficient) >= minRegressionCoefficient.get();
                    }).forEachOrdered(line -> {
                        try {
                            writer.write(line);
                            writer.newLine();
                        } catch (IOException e) {
                            throw new UncheckedIOException(e);
                        }
                    });
                }
                return true;
            }));
        }};
    }
}
