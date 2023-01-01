package org.exbio.tfprio.steps.plots;

import org.apache.commons.lang3.tuple.Pair;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;

public class GroupStages extends ExecutableStep {
    public final Map<String, Map<Double, Pair<OutputFile, OutputFile>>> outputFiles = new HashMap<>();
    private final RequiredConfig<Map<String, String>> sameStages =
            new RequiredConfig<>(Configs.inputConfigs.sameStages);
    private final RequiredConfig<List<Double>> thresholds = new RequiredConfig<>(Configs.plots.thresholds);
    private final Map<OutputFile, Collection<InputFile>> bridge = new HashMap<>();
    private final Map<InputFile, String> inputFilePairing = new HashMap<>();

    public GroupStages(Map<String, Map<String, Map<Double, OutputFile>>> coefficientFiles) {
        super(false, coefficientFiles.values().stream().flatMap(x -> x.values().stream()).flatMap(
                col -> col.values().stream()).toList());

        coefficientFiles.values().stream().map(Map::keySet).flatMap(Collection::stream).distinct().forEach(hm -> {
            OutputFile hmDir = new OutputFile(outputDirectory, hm);

            thresholds.get().forEach(threshold -> {
                OutputFile thresholdDir = new OutputFile(hmDir, threshold.toString());
                outputFiles.computeIfAbsent(hm, s -> new HashMap<>()).put(threshold,
                        Pair.of(addOutput(thresholdDir, "same.tsv"), addOutput(thresholdDir, "different.tsv")));
            });
        });


        coefficientFiles.forEach((pairing, hmMap) -> {
            String[] pairingSplit = pairing.split("_");
            String group1 = pairingSplit[0];
            String group2 = pairingSplit[1];

            boolean same = sameStages.get().get(group1).equals(sameStages.get().get(group2));

            OutputFile inPairing = new OutputFile(inputDirectory, pairing);

            hmMap.forEach((hm, thresholdMap) -> {
                OutputFile hmDir = new OutputFile(inPairing, hm);
                thresholdMap.forEach((threshold, input) -> {
                    InputFile inThreshold = addInput(hmDir, input);
                    OutputFile outputFile = same ? outputFiles.get(hm).get(threshold).getLeft() :
                            outputFiles.get(hm).get(threshold).getRight();
                    bridge.computeIfAbsent(outputFile, s -> new HashSet<>()).add(inThreshold);
                    inputFilePairing.put(inThreshold, pairing);
                });
            });
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((outputFile, inputFiles) -> add(() -> {
                Map<String, Map<String, Double>> tfPairingValues = new HashMap<>();

                inputFiles.forEach(inputFile -> {
                    String pairing = inputFilePairing.get(inputFile);

                    try (BufferedReader reader = new BufferedReader(new FileReader(inputFile))) {
                        reader.lines().skip(1).map(line -> line.split("\t")).map(
                                split -> Pair.of(split[0], Double.parseDouble(split[1]))).forEach(pair -> {
                            tfPairingValues.computeIfAbsent(pair.getKey(), s -> new HashMap<>()).put(pairing,
                                    pair.getValue());
                        });
                    } catch (IOException e) {
                        logger.error("Error reading file", e);
                    }
                });

                List<String> pairings = tfPairingValues.values().stream().map(Map::keySet).flatMap(
                        Collection::stream).distinct().sorted().toList();

                try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
                    writer.write("TF\t" + String.join("\t", pairings));
                    writer.newLine();

                    tfPairingValues.forEach((tf, pairingValues) -> {
                        try {
                            writer.write(tf + "\t" + pairings.stream().map(
                                    pairing -> pairingValues.getOrDefault(pairing, Double.NaN)).map(
                                    Object::toString).collect(java.util.stream.Collectors.joining("\t")));
                            writer.newLine();
                        } catch (IOException e) {
                            logger.error("Error writing file", e);
                        }
                    });
                } catch (IOException e) {
                    logger.error("Error writing file", e);
                }

                return true;
            }));
        }};
    }
}
