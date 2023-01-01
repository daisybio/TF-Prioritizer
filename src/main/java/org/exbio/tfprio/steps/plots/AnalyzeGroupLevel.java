package org.exbio.tfprio.steps.plots;

import org.apache.commons.lang3.tuple.Pair;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.OptionalConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static java.util.stream.Collectors.toMap;
import static org.exbio.pipejar.util.FileManagement.readLines;

public class AnalyzeGroupLevel extends ExecutableStep {
    private final InputFile ensgSymbolFile;
    private final InputFile analyzableTfs;
    private final Map<String, InputFile> groupTPMFiles;
    private final Map<String, InputFile> groupMeanCountFiles;

    private final Map<InputFile, OutputFile> bridge = new HashMap<>();
    private final OptionalConfig<Integer> cutoffGroups = new OptionalConfig<>(Configs.plots.cutoffGroups, false);
    private final OptionalConfig<Double> cutoffGroupCounts =
            new OptionalConfig<>(Configs.plots.cutoffGroupCounts, false);
    private final OptionalConfig<Double> cutoffTPM = new OptionalConfig<>(Configs.plots.cutoffTPM, false);


    public AnalyzeGroupLevel(Map<String, OutputFile> groupMeanTPM, Map<String, OutputFile> groupMeanCount,
                             Map<String, Map<Double, Pair<OutputFile, OutputFile>>> hmThresholdGroupedStages,
                             OutputFile ensgSymbolFile, OutputFile analyzableTfs) {
        super(false, groupMeanTPM.values(), groupMeanCount.values(),
                hmThresholdGroupedStages.values().stream().flatMap(x -> x.values().stream()).flatMap(
                        pair -> Stream.of(pair.getLeft(), pair.getRight())).toList(),
                List.of(ensgSymbolFile, analyzableTfs));

        this.ensgSymbolFile = addInput(ensgSymbolFile);
        this.analyzableTfs = addInput(analyzableTfs);

        OutputFile inputTpm = new OutputFile(inputDirectory, "tpm");
        OutputFile inputCount = new OutputFile(inputDirectory, "count");
        OutputFile inputGroupedStages = new OutputFile(inputDirectory, "groupedStages");

        groupTPMFiles = groupMeanTPM.entrySet().stream().collect(
                toMap(Map.Entry::getKey, e -> addInput(inputTpm, e.getValue())));
        groupMeanCountFiles = groupMeanCount.entrySet().stream().collect(
                toMap(Map.Entry::getKey, e -> addInput(inputCount, e.getValue())));

        hmThresholdGroupedStages.forEach((hm, thresholdMap) -> {
            OutputFile inHm = new OutputFile(inputGroupedStages, hm);
            OutputFile outHm = new OutputFile(outputDirectory, hm);
            thresholdMap.forEach((threshold, pair) -> {
                OutputFile inThreshold = new OutputFile(inHm, threshold.toString());
                OutputFile outThreshold = new OutputFile(outHm, threshold.toString());

                InputFile inputSame = addInput(inThreshold, pair.getLeft());
                InputFile inputDifferent = addInput(inThreshold, pair.getRight());

                OutputFile outputSame = addOutput(outThreshold, "same.tsv");
                OutputFile outputDifferent = addOutput(outThreshold, "different.tsv");

                bridge.put(inputSame, outputSame);
                bridge.put(inputDifferent, outputDifferent);
            });
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            final Map<String, String> ensgSymbol;
            try {
                ensgSymbol = readLines(ensgSymbolFile).stream().map(line -> line.split("\t")).collect(
                        toMap(line -> line[0], line -> line.length > 1 ? line[1].toUpperCase() : line[0]));
            } catch (IOException e) {
                throw new RuntimeException(e);
            }

            Map<String, Map<String, Double>> groupGeneNameTpm = readExpression(groupTPMFiles, ensgSymbol);
            Map<String, Map<String, Double>> groupGeneNameCount = readExpression(groupMeanCountFiles, ensgSymbol);

            final List<String> analyzableTfIds;
            try {
                analyzableTfIds = readLines(analyzableTfs);
            } catch (IOException e) {
                throw new RuntimeException(e);
            }

            List<String> groupOrder = groupGeneNameCount.keySet().stream().sorted().toList();

            Map<String, Map<String, Double>> groupTfCount = groupGeneNameCount.entrySet().stream().collect(
                    toMap(Map.Entry::getKey, geneNameCount -> analyzableTfIds.stream().map(String::toUpperCase).collect(
                            toMap(Function.identity(), tfId -> tfId.contains("::") ?
                                    Arrays.stream(tfId.split("::")).mapToDouble(
                                            tfPart -> geneNameCount.getValue().getOrDefault(tfPart, 0.0)).sum() :
                                    geneNameCount.getValue().getOrDefault(tfId, 0.0)))));

            bridge.forEach((input, output) -> add(() -> {
                try (BufferedReader reader = new BufferedReader(new FileReader(input));
                     BufferedWriter writer = new BufferedWriter(new FileWriter(output))) {
                    writer.write("TF\t" + String.join("\t", groupOrder));
                    writer.newLine();

                    reader.lines().skip(1).map(line -> line.split("\t"))
                          // Split into tf und coefficient values
                          .map(split -> Pair.of(split[0],
                                  // Empty strings mean that a tf was not found in a group
                                  Stream.of(split).skip(1).filter(s -> !s.isEmpty()).toList()))
                          // Keep only entries with more than cutoffGroups entries, if cutoffGroups is set
                          .filter(pair -> !cutoffGroups.isSet() || pair.getRight().size() > cutoffGroups.get())
                          // Remove coefficients
                          .map(Pair::getLeft).map(String::toUpperCase)
                          // Map tf to a Pair of tf and group-count map
                          .forEach(tf -> {
                              // Get the count for the tf in each group
                              Map<String, Double> groupCounts = groupTfCount.entrySet().stream().collect(
                                      toMap(Map.Entry::getKey, entry -> entry.getValue().getOrDefault(tf, 0.0)));

                              Double countSum = groupCounts.values().stream().mapToDouble(Double::doubleValue).sum();
                              Double tpmSum = groupGeneNameTpm.values().stream().mapToDouble(
                                      stringDoubleMap -> stringDoubleMap.getOrDefault(tf, 0.0)).sum();

                              if (cutoffGroupCounts.isSet() && countSum < cutoffGroupCounts.get()) {
                                  return;
                              }
                              if (cutoffTPM.isSet() && tpmSum < cutoffTPM.get()) {
                                  return;
                              }

                              if (countSum <= 0) {
                                  return;
                              }

                              try {
                                  writer.write(tf + "\t" + String.join("\t", groupOrder.stream().map(
                                          group -> groupCounts.getOrDefault(group, 0.0).toString()).toList()));
                                  writer.newLine();
                              } catch (IOException e) {
                                  throw new RuntimeException(e);
                              }
                          });
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
                return true;
            }));
        }};
    }

    private Map<String, Map<String, Double>> readExpression(Map<String, InputFile> expressionFiles,
                                                            Map<String, String> ensgSymbol) {
        return expressionFiles.entrySet().stream().collect(toMap(Map.Entry::getKey, e -> {
            try {
                return readLines(e.getValue()).stream().skip(1).map(line -> line.split("\t")).map(
                        split -> Pair.of(ensgSymbol.getOrDefault(split[0], split[0]),
                                Double.parseDouble(split[1]))).collect(
                        Collectors.groupingBy(Pair::getLeft, Collectors.summingDouble(Pair::getValue)));
            } catch (IOException ex) {
                throw new RuntimeException(ex);
            }
        }));
    }
}
