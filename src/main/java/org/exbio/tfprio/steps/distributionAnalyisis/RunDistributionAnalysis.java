package org.exbio.tfprio.steps.distributionAnalyisis;

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
import java.util.stream.IntStream;

import static org.exbio.pipejar.util.FileManagement.readLines;

public class RunDistributionAnalysis extends ExecutableStep {
    public final Map<String, Map<String, OutputFile>> outputFiles = new HashMap<>();
    private final InputFile preprocessed;
    private final InputFile ensgSymbol;
    private final Map<String, InputFile> groupMeanCounts = new HashMap<>();
    private final Map<String, InputFile> pairingLog2fc = new HashMap<>();
    private final Map<String, Map<String, InputFile>> pairingHmRegressionCoefficients = new HashMap<>();
    private final Map<String, Map<String, InputFile>> groupHmAffinityFiles = new HashMap<>();
    private final OptionalConfig<Boolean> scoreIncludeCounts =
            new OptionalConfig<>(Configs.distributionAnalysis.scoreIncludeCounts, false);

    public RunDistributionAnalysis(OutputFile preprocessed, OutputFile ensgSymbol,
                                   Map<String, OutputFile> groupMeanCounts, Map<String, OutputFile> pairingLog2fc,
                                   Map<String, Map<String, OutputFile>> pairingHmRegressionCoefficients,
                                   Map<String, Map<String, OutputFile>> groupHmAffinityFiles) {
        super(false, pairingLog2fc.values(), groupMeanCounts.values(),
                pairingHmRegressionCoefficients.values().stream().flatMap(entry -> entry.values().stream()).toList(),
                groupHmAffinityFiles.values().stream().flatMap(entry -> entry.values().stream()).toList(),
                List.of(preprocessed, ensgSymbol));

        this.preprocessed = addInput(preprocessed);
        this.ensgSymbol = addInput(ensgSymbol);

        OutputFile meanCountsDir = new OutputFile(inputDirectory, "counts");
        groupMeanCounts.forEach(
                (group, meanCounts) -> this.groupMeanCounts.put(group, addInput(meanCountsDir, meanCounts)));

        OutputFile diffExpression = new OutputFile(inputDirectory, "diffExpression");
        pairingLog2fc.forEach((pairing, log2fc) -> this.pairingLog2fc.put(pairing, addInput(diffExpression, log2fc)));

        OutputFile regressionCoefficients = new OutputFile(inputDirectory, "regressionCoefficients");
        pairingHmRegressionCoefficients.forEach((pairing, hmMap) -> {
            OutputFile pairingDir = new OutputFile(regressionCoefficients, pairing);
            OutputFile outPairing = new OutputFile(outputDirectory, pairing);
            hmMap.forEach((hm, regressionCoefficientFile) -> {
                this.pairingHmRegressionCoefficients.computeIfAbsent(pairing, s -> new HashMap<>()).put(hm,
                        addInput(pairingDir, regressionCoefficientFile));
                this.outputFiles.computeIfAbsent(pairing, s -> new HashMap<>()).put(hm, addOutput(outPairing, hm));
            });
        });

        OutputFile affinityFiles = new OutputFile(inputDirectory, "affinityFiles");
        groupHmAffinityFiles.forEach((group, hmMap) -> {
            OutputFile groupDir = new OutputFile(affinityFiles, group);

            hmMap.forEach(
                    (hm, affinityFile) -> this.groupHmAffinityFiles.computeIfAbsent(group, s -> new HashMap<>()).put(hm,
                            addInput(groupDir, affinityFile)));
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            Map<String, Collection<String>> tfHms;
            try {
                tfHms = readLines(preprocessed).stream().skip(1).map(line -> line.split("\t")).map(
                        split -> Pair.of(split[0], split[1])).collect(Collectors.groupingBy(Pair::getKey,
                        Collectors.mapping(Pair::getValue, Collectors.toCollection(HashSet::new))));
            } catch (IOException e) {
                throw new RuntimeException(e);
            }


            Map<String, Map<String, Double>> groupGeneCount = readStateEntryTable(groupMeanCounts);
            Map<String, Map<String, Double>> pairingGeneLog2fc = readStateEntryTable(pairingLog2fc);
            Map<String, Double> pairingAverageLog2fc = pairingGeneLog2fc.entrySet().stream().collect(
                    Collectors.toMap(Map.Entry::getKey, entry -> entry.getValue().values().stream().mapToDouble(
                            Double::doubleValue).average().orElse(0.0)));

            Map<String, Map<String, Map<String, Double>>> pairingHmTfRegressionCoefficients =
                    pairingHmRegressionCoefficients.entrySet().stream().map(pairingEntry -> {
                        String pairing = pairingEntry.getKey();
                        Map<String, InputFile> hmRegressionCoefficientFile = pairingEntry.getValue();

                        Map<String, Map<String, Double>> hmTfRegressionCoefficient =
                                readStateEntryTable(hmRegressionCoefficientFile);
                        return Pair.of(pairing, hmTfRegressionCoefficient);
                    }).collect(Collectors.toMap(Pair::getKey, Pair::getValue));

            Set<String> tfs = tfHms.keySet();

            Map<String, Map<String, Map<String, Map<String, Double>>>> groupHmTfGeneAffinity =
                    groupHmAffinityFiles.entrySet().stream().map(groupEntry -> {
                        String group = groupEntry.getKey();

                        Map<String, InputFile> hmAffinityFile = groupEntry.getValue();

                        Map<String, Map<String, Map<String, Double>>> hmTfGeneAffinity =
                                hmAffinityFile.entrySet().stream().map(hmEntry -> {
                                    String hm = hmEntry.getKey();
                                    InputFile affinityFile = hmEntry.getValue();

                                    Map<String, Map<String, Double>> tfGeneAffinity;

                                    try (BufferedReader reader = new BufferedReader(new FileReader(affinityFile))) {
                                        String[] header = reader.readLine().split("\t");

                                        Map<String, Integer> tfIndex = IntStream.range(0, header.length).filter(
                                                i -> tfs.contains(header[i])).boxed().collect(
                                                Collectors.toMap(i -> header[i], Function.identity()));

                                        Map<String, Map<String, Double>> geneTfAffinity =
                                                reader.lines().map(line -> line.split("\t")).map(split -> {
                                                    String gene = split[0];
                                                    Map<String, Double> tfAffinity = tfIndex.entrySet().stream().map(
                                                            tfEntry -> Pair.of(tfEntry.getKey(), Double.parseDouble(
                                                                    split[tfEntry.getValue()]))).collect(
                                                            Collectors.toMap(Pair::getKey, Pair::getValue));
                                                    return Pair.of(gene, tfAffinity);
                                                }).collect(Collectors.toMap(Pair::getKey, Pair::getValue));

                                        tfGeneAffinity = tfIndex.keySet().stream().collect(
                                                Collectors.toMap(Function.identity(),
                                                        tf -> geneTfAffinity.entrySet().stream().filter(
                                                                geneEntry -> geneEntry.getValue().containsKey(tf)).map(
                                                                geneEntry -> Pair.of(geneEntry.getKey(),
                                                                        geneEntry.getValue().get(tf))).collect(
                                                                Collectors.toMap(Pair::getKey, Pair::getValue))));
                                    } catch (IOException e) {
                                        throw new RuntimeException(e);
                                    }


                                    return Pair.of(hm, tfGeneAffinity);
                                }).collect(Collectors.toMap(Pair::getKey, Pair::getValue));
                        return Pair.of(group, hmTfGeneAffinity);
                    }).collect(Collectors.toMap(Pair::getKey, Pair::getValue));

            tfHms.forEach((tf, hms) -> add(() -> {
                pairingAverageLog2fc.forEach((pairing, averageLog2fc) -> {
                    String[] split = pairing.split("_");
                    String group1 = split[0];
                    String group2 = split[1];

                    hms.forEach(hm -> {
                        // Here the available genes were extracted from tgene files (line 231-250)
                        // Not implemented, because Tgene files are not available
                        // TODO implement

                        double regressionCoefficient =
                                pairingHmTfRegressionCoefficients.getOrDefault(pairing, new HashMap<>()).getOrDefault(
                                        hm, new HashMap<>()).getOrDefault(tf, -1.);

                        if (regressionCoefficient <= 0.0) {
                            return;
                        }

                        Map<String, Double> targetGeneAffinity1 =
                                groupHmTfGeneAffinity.getOrDefault(group1, new HashMap<>()).getOrDefault(hm,
                                        new HashMap<>()).getOrDefault(tf, new HashMap<>());
                        Map<String, Double> targetGeneAffinity2 =
                                groupHmTfGeneAffinity.getOrDefault(group2, new HashMap<>()).getOrDefault(hm,
                                        new HashMap<>()).getOrDefault(tf, new HashMap<>());

                        Set<String> targetGenes = new HashSet<>(targetGeneAffinity1.keySet());
                        targetGenes.retainAll(targetGeneAffinity2.keySet());

                        OutputFile outDir = outputFiles.get(pairing).get(hm);
                        File outFile = new File(outDir, tf + ".tsv");

                        if (targetGenes.isEmpty()) {
                            return;
                        }

                        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outFile))) {
                            writer.write("Gene\tTfTgScore\tHM\tPairing\tTF\tRegressionCoefficient");
                            writer.newLine();

                            targetGenes.forEach(targetGene -> {
                                double affinity1 = targetGeneAffinity1.get(targetGene);
                                double affinity2 = targetGeneAffinity2.get(targetGene);

                                double geneCount1 =
                                        groupGeneCount.getOrDefault(group1, new HashMap<>()).getOrDefault(targetGene,
                                                0.0);
                                double geneCount2 =
                                        groupGeneCount.getOrDefault(group2, new HashMap<>()).getOrDefault(targetGene,
                                                0.0);


                                double tgScore1 = averageLog2fc * affinity1;
                                double tgScore2 = averageLog2fc * affinity2;

                                if (scoreIncludeCounts.isSet() && scoreIncludeCounts.get()) {
                                    tgScore1 *= geneCount1;
                                    tgScore2 *= geneCount2;
                                }

                                double tgScore = Math.abs(tgScore1) + Math.abs(tgScore2);
                                double tfTgScore = Math.abs(tgScore * regressionCoefficient);

                                try {
                                    writer.write(
                                            targetGene + "\t" + tfTgScore + "\t" + hm + "\t" + pairing + "\t" + tf +
                                                    "\t" + regressionCoefficient);
                                    writer.newLine();
                                } catch (IOException e) {
                                    throw new RuntimeException(e);
                                }
                            });
                        } catch (IOException e) {
                            throw new RuntimeException(e);
                        }
                    });
                });
                return true;
            }));
        }};
    }

    private Map<String, Map<String, Double>> readStateEntryTable(Map<String, InputFile> stateFileMap) {
        return stateFileMap.entrySet().stream().map(entry -> {
            String state = entry.getKey();
            InputFile stateFile = entry.getValue();

            Map<String, Double> geneValue;
            try {
                geneValue = readLines(stateFile).stream().skip(1).map(line -> line.split("\t")).collect(
                        Collectors.toMap(split -> split[0], split -> Double.parseDouble(split[1])));
            } catch (IOException e) {
                throw new RuntimeException(e);
            }

            return Pair.of(state, geneValue);
        }).collect(Collectors.toMap(Pair::getKey, Pair::getValue));
    }
}
