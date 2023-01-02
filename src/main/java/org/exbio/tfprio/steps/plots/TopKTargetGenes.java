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
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static org.exbio.pipejar.util.FileManagement.makeSureFileExists;

public class TopKTargetGenes extends ExecutableStep {
    public final Map<String, Map<Double, Pair<OutputFile, OutputFile>>> outputFiles = new HashMap<>();
    private final Map<OutputFile, Pair<InputFile, Map<String, InputFile>>> bridge = new HashMap<>();
    private final RequiredConfig<Integer> topKTargetGenes = new RequiredConfig<>(Configs.plots.topKTargetGenes);

    public TopKTargetGenes(Map<String, Map<Double, Pair<OutputFile, OutputFile>>> hmThresholdGrouped,
                           Map<String, Map<String, OutputFile>> groupHmMeanAffinities) {
        super(false, hmThresholdGrouped.values().stream().flatMap(x -> x.values().stream()).flatMap(
                        pair -> Stream.of(pair.getLeft(), pair.getRight())).toList(),
                groupHmMeanAffinities.values().stream().flatMap(sub -> sub.values().stream()).toList());

        Set<String> groups = hmThresholdGrouped.keySet();

        OutputFile inCoefficients = new OutputFile(inputDirectory, "coefficients");
        OutputFile inAffinities = new OutputFile(inputDirectory, "meanAffinities");

        hmThresholdGrouped.forEach((hm, thresholdMap) -> {
            OutputFile inCoefficientsHm = new OutputFile(inCoefficients, hm);
            OutputFile inAffinitiesHm = new OutputFile(inAffinities, hm);
            OutputFile outHm = new OutputFile(outputDirectory, hm);

            Map<String, InputFile> groupMeanAffinities = groupHmMeanAffinities.entrySet().stream().map(
                    entry -> Pair.of(entry.getKey(), entry.getValue().get(hm))).filter(
                    pair -> !Objects.isNull(pair.getValue())).map(
                    pair -> Pair.of(pair.getKey(), addInput(inAffinitiesHm, pair.getValue()))).collect(
                    Collectors.toMap(Pair::getKey, Pair::getValue));


            thresholdMap.forEach((threshold, pair) -> {
                OutputFile inThreshold = new OutputFile(inCoefficientsHm, threshold.toString());
                OutputFile outThreshold = new OutputFile(outHm, threshold.toString());

                InputFile inSame = addInput(inThreshold, pair.getLeft());
                InputFile inDifferent = addInput(inThreshold, pair.getRight());

                OutputFile outSame = addOutput(outThreshold, "same");
                OutputFile outDifferent = addOutput(outThreshold, "different");

                bridge.put(outSame, Pair.of(inSame, groupMeanAffinities));
                bridge.put(outDifferent, Pair.of(inDifferent, groupMeanAffinities));
            });
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((output, pair) -> add(() -> {
                InputFile inputCoefficients = pair.getLeft();
                Map<String, InputFile> groupMeanAffinities = pair.getRight();

                Map<String, Collection<String>> groupTfs;

                try (var reader = new BufferedReader(new FileReader(inputCoefficients))) {
                    List<String> pairings = Arrays.stream(reader.readLine().split("\t")).skip(1).toList();

                    Map<String, Collection<String>> tfPairings = reader.lines().map(line -> line.split("\t")).map(
                            split -> Pair.of(split[0], List.of(split).subList(1, split.length))).map(
                            p -> Pair.of(p.getLeft(), IntStream.range(0, p.getRight().size()).filter(
                                    i -> !p.getRight().get(i).isEmpty()).mapToObj(pairings::get).toList())).collect(
                            HashMap::new, (map, p) -> map.put(p.getLeft(), p.getRight()), HashMap::putAll);

                    Map<String, Collection<String>> tfGroups = tfPairings.entrySet().stream().map(
                            entry -> Pair.of(entry.getKey(),
                                    entry.getValue().stream().map(pairing -> pairing.split("_")).flatMap(
                                            Arrays::stream).distinct().toList())).collect(
                            Collectors.toMap(Pair::getKey, Pair::getValue));

                    groupTfs = tfGroups.entrySet().stream().flatMap(
                            entry -> entry.getValue().stream().map(group -> Pair.of(group, entry.getKey()))).collect(
                            Collectors.groupingBy(Pair::getKey,
                                    Collectors.mapping(Pair::getValue, Collectors.toCollection(HashSet::new))));
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }

                groupTfs.forEach((group, tfs) -> {
                    File outputGroup = new File(output, group);
                    InputFile inputAffinities = groupMeanAffinities.get(group);

                    try (var reader = new BufferedReader(new FileReader(inputAffinities))) {
                        String[] header = reader.readLine().split("\t");
                        Map<String, Integer> tfIndices = IntStream.range(1, header.length).mapToObj(
                                i -> Pair.of(header[i].toUpperCase(), i)).filter(p -> tfs.contains(p.getKey())).collect(
                                Collectors.toMap(Pair::getKey, Pair::getValue));

                        Map<String, Collection<Pair<String, Double>>> tfGeneAffinity =
                                reader.lines().map(line -> line.split("\t")).flatMap(split -> {
                                    String gene = split[0];
                                    return tfIndices.entrySet().stream().map(entry -> {
                                        double affinity = Double.parseDouble(split[entry.getValue()]);
                                        return Pair.of(entry.getKey(), Pair.of(gene, affinity));
                                    });
                                }).collect(Collectors.groupingBy(Pair::getKey,
                                        Collectors.mapping(Pair::getValue, Collectors.toCollection(HashSet::new))));

                        Map<String, List<Pair<String, Double>>> tfTopGeneAffinity =
                                tfGeneAffinity.entrySet().stream().map(entry -> Pair.of(entry.getKey(),
                                        entry.getValue().stream().sorted(
                                                (a, b) -> b.getValue().compareTo(a.getValue())).limit(
                                                topKTargetGenes.get()).toList())).collect(
                                        Collectors.toMap(Pair::getKey, Pair::getValue));

                        tfTopGeneAffinity.forEach((tf, topGeneAffinity) -> {
                            File outputTf = new File(outputGroup, tf + ".tsv");
                            try {
                                makeSureFileExists(outputTf);
                                try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputTf))) {
                                    writer.write("Gene\tAffinity");
                                    writer.newLine();

                                    topGeneAffinity.forEach(p -> {
                                        try {
                                            writer.write(p.getKey() + "\t" + p.getValue());
                                            writer.newLine();
                                        } catch (IOException e) {
                                            throw new RuntimeException(e);
                                        }
                                    });
                                }
                            } catch (IOException e) {
                                throw new RuntimeException(e);
                            }
                        });
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                });

                return true;
            }));
        }};
    }
}
