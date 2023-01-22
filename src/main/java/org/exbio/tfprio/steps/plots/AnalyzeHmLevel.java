package org.exbio.tfprio.steps.plots;

import org.apache.commons.lang3.tuple.Pair;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.OptionalConfig;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class AnalyzeHmLevel extends ExecutableStep<Configs> {
    public final Map<Double, Pair<OutputFile, OutputFile>> outputFiles = new HashMap<>();
    private final Map<OutputFile, Collection<InputFile>> bridge = new HashMap<>();
    private final RequiredConfig<List<Double>> thresholds = new RequiredConfig<>(configs.plots.thresholds);

    private final OptionalConfig<Integer> cutoffHms = new OptionalConfig<>(configs.plots.cutoffHms, false);

    public AnalyzeHmLevel(Configs configs,
                          Map<String, Map<Double, Pair<OutputFile, OutputFile>>> hmThresholdGroupLevel) {
        super(configs, false, hmThresholdGroupLevel.values().stream().flatMap(x -> x.values().stream()).flatMap(
                pair -> Stream.of(pair.getLeft(), pair.getRight())).toList());

        thresholds.get().forEach(threshold -> {
            OutputFile thresholdDir = new OutputFile(outputDirectory, threshold.toString());
            outputFiles.put(threshold,
                    Pair.of(addOutput(thresholdDir, "same.tsv"), addOutput(thresholdDir, "different.tsv")));
        });

        hmThresholdGroupLevel.forEach((hm, thresholdMap) -> {
            OutputFile inHm = new OutputFile(inputDirectory, hm);
            thresholdMap.forEach((threshold, pair) -> {
                OutputFile inThreshold = new OutputFile(inHm, threshold.toString());

                InputFile inputSame = addInput(inThreshold, pair.getLeft());
                InputFile inputDiff = addInput(inThreshold, pair.getRight());

                bridge.computeIfAbsent(outputFiles.get(threshold).getLeft(), s -> new HashSet<>()).add(inputSame);
                bridge.computeIfAbsent(outputFiles.get(threshold).getRight(), s -> new HashSet<>()).add(inputDiff);
            });
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((outputFile, inputFiles) -> add(() -> {
                Set<String> groups = new HashSet<>();

                Collection<Map<String, Map<String, Double>>> tfGroupCounts = inputFiles.stream().map(inputFile -> {
                    try (BufferedReader reader = new BufferedReader(new FileReader(inputFile))) {
                        String[] headerSplit = reader.readLine().split("\t");
                        groups.addAll(Arrays.stream(headerSplit).skip(1).collect(Collectors.toSet()));

                        Map<String, Map<String, Double>> tfGroupCount =
                                reader.lines().map(line -> line.split("\t")).map(
                                        split -> Pair.of(split[0], List.of(split).subList(1, split.length))).map(
                                        pair -> {
                                            Map<String, Double> groupCount =
                                                    pair.getRight().stream().map(Double::parseDouble).collect(
                                                            HashMap::new,
                                                            (map, count) -> map.put(headerSplit[map.size() + 1], count),
                                                            HashMap::putAll);
                                            return Pair.of(pair.getLeft(), groupCount);
                                        }).collect(HashMap::new,
                                        (map, pair) -> map.put(pair.getLeft(), pair.getRight()), HashMap::putAll);
                        return tfGroupCount;
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                }).collect(Collectors.toSet());

                Collection<String> tfs = tfGroupCounts.stream().flatMap(map -> map.keySet().stream()).filter(
                        tf -> !cutoffHms.isSet() || tfGroupCounts.stream().filter(map -> map.containsKey(tf)).count() >
                                cutoffHms.get()).collect(Collectors.toSet());

                List<String> groupOrder = groups.stream().sorted().toList();

                try (var writer = new BufferedWriter(new FileWriter(outputFile))) {
                    writer.write("TF\t" + String.join("\t", groupOrder));
                    writer.newLine();
                    tfs.forEach(tf -> {
                        try {
                            writer.write(tf + "\t" + String.join("\t", groupOrder.stream().map(
                                    group -> tfGroupCounts.stream().filter(map -> map.containsKey(tf)).mapToDouble(
                                            map -> map.get(tf).get(group)).sum()).map(String::valueOf).toList()));
                            writer.newLine();
                        } catch (IOException e) {
                            throw new RuntimeException(e);
                        }
                    });
                }

                return true;
            }));
        }};
    }
}
