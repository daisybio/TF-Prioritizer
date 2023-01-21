package org.exbio.tfprio.steps.Dynamite;

import org.apache.commons.lang3.tuple.Pair;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;

public class FilterRegressionCoefficients extends ExecutableStep<Configs> {
    public final Map<String, Map<String, Map<Double, OutputFile>>> outputFiles = new HashMap<>();
    private final Map<InputFile, Map<Double, OutputFile>> bridge = new HashMap<>();

    private final RequiredConfig<List<Double>> thresholds = new RequiredConfig<>(configs.plots.thresholds);

    public FilterRegressionCoefficients(Configs configs, Map<String, Map<String, OutputFile>> regressionCoefficients) {
        super(configs, false, regressionCoefficients.values().stream().flatMap(
                stringCollectionMap -> stringCollectionMap.values().stream()).toList());

        regressionCoefficients.forEach((pairing, hmMap) -> {
            OutputFile inputPairing = new OutputFile(inputDirectory, pairing);
            OutputFile outputPairing = new OutputFile(outputDirectory, pairing);

            hmMap.forEach((hm, inDir) -> {
                InputFile inputDir = addInput(inputPairing, inDir);
                OutputFile outputHm = new OutputFile(outputPairing, hm);
                thresholds.get().forEach(threshold -> {
                    OutputFile outputThreshold = addOutput(outputHm, threshold + ".tsv");

                    bridge.computeIfAbsent(inputDir, s -> new HashMap<>()).put(threshold, outputThreshold);
                    outputFiles.computeIfAbsent(pairing, s -> new HashMap<>()).computeIfAbsent(hm,
                            s -> new HashMap<>()).put(threshold, outputThreshold);
                });
            });
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((extractedCoefficients, hmMap) -> hmMap.forEach((threshold, outputThreshold) -> add(() -> {
                try (BufferedReader reader = new BufferedReader(new FileReader(extractedCoefficients));
                     BufferedWriter writer = new BufferedWriter(new FileWriter(outputThreshold))) {

                    writer.write(reader.readLine());
                    writer.newLine();

                    reader.lines().map(line -> line.split("\t")).map(
                            split -> Pair.of(split[0], Double.parseDouble(split[1]))).filter(pair -> {
                        double value = pair.getValue();
                        return value >= threshold || value <= -threshold;
                    }).sorted(Comparator.comparingDouble(Pair::getValue)).forEach(pair -> {
                        try {
                            writer.write(pair.getKey() + "\t" + pair.getValue());
                            writer.newLine();
                        } catch (IOException e) {
                            throw new RuntimeException(e);
                        }
                    });
                }

                return true;
            })));
        }};
    }
}
