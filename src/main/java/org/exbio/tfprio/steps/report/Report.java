package org.exbio.tfprio.steps.report;

import org.apache.commons.lang3.tuple.Pair;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.json.JSONObject;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import static java.util.stream.Collectors.toMap;
import static org.exbio.pipejar.util.FileManagement.readLines;

public class Report extends ExecutableStep {
    public final OutputFile outputFile = addOutput("report");
    private final InputFile ensgSymbol;
    private final Map<String, InputFile> hmRanks = new HashMap<>();
    private final Map<String, InputFile> pairingDeseq = new HashMap<>();
    private final Map<String, InputFile> groupMeanExpression = new HashMap<>();
    private final Map<String, InputFile> groupTpm = new HashMap<>();

    public Report(OutputFile ensgSymbol, Map<String, OutputFile> hmRanks, Map<String, OutputFile> deseqResults,
                  Map<String, OutputFile> groupMeanExpression, Map<String, OutputFile> groupTpm) {
        super(false, hmRanks.values(), deseqResults.values(), groupMeanExpression.values(), groupTpm.values(),
                List.of(ensgSymbol));

        this.ensgSymbol = addInput(ensgSymbol);
        OutputFile rankDir = new OutputFile(inputDirectory, "ranks");
        hmRanks.forEach((hm, rankFile) -> this.hmRanks.put(hm, addInput(rankDir, rankFile)));

        OutputFile deseqDir = new OutputFile(inputDirectory, "deseq");
        deseqResults.forEach((pairing, deseqFile) -> this.pairingDeseq.put(pairing, addInput(deseqDir, deseqFile)));

        OutputFile meanDir = new OutputFile(inputDirectory, "meanExpression");
        groupMeanExpression.forEach(
                (group, meanFile) -> this.groupMeanExpression.put(group, addInput(meanDir, meanFile)));

        OutputFile tpmDir = new OutputFile(inputDirectory, "tpm");
        groupTpm.forEach((group, tpmFile) -> this.groupTpm.put(group, addInput(tpmDir, tpmFile)));
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                Map<String, List<String>> symbolEnsgMap =
                        readLines(ensgSymbol).stream().map(line -> line.split("\t")).filter(
                                split -> split.length > 1).collect(Collectors.groupingBy(split -> split[1],
                                Collectors.mapping(split -> split[0], Collectors.toList())));

                logger.trace("Fetching tf symbols");

                Set<String> tfSymbols = hmRanks.values().stream().flatMap(rankFile -> {
                    try {
                        return readLines(rankFile).stream().skip(1);
                    } catch (IOException e) {
                        throw new UncheckedIOException(e);
                    }
                }).map(line -> line.split("\t")[0]).collect(Collectors.toSet());

                logger.trace("Creating tf groups");

                Collection<TfGroup> groups = tfSymbols.stream().map(tfSymbol -> {
                    String[] tfParts = tfSymbol.split("::");

                    Collection<TranscriptionFactor> transcriptionFactors = Arrays.stream(tfParts).map(tfPart -> {
                        List<String> ensgs = symbolEnsgMap.get(tfPart);
                        return new TranscriptionFactor(tfPart, ensgs);
                    }).toList();

                    return new TfGroup(tfSymbol, transcriptionFactors);
                }).collect(Collectors.toSet());

                logger.trace("Reading log2 fold changes");

                pairingDeseq.forEach((pairing, deseqFile) -> {
                    Map<String, Double> ensgLog2fc;
                    try {
                        ensgLog2fc = readLines(deseqFile).stream().skip(1).map(line -> line.split("\t")).map(
                                split -> Pair.of(split[0].replace("\"", ""), Double.parseDouble(split[1]))).collect(
                                toMap(Pair::getKey, Pair::getValue));
                    } catch (IOException e) {
                        throw new UncheckedIOException(e);
                    }
                    groups.forEach(group -> group.getTranscriptionFactors().forEach(tf -> tf.setLog2fc(pairing,
                            tf.getEnsgs().stream().filter(ensgLog2fc::containsKey).mapToDouble(
                                    ensgLog2fc::get).average().orElse(0.0))));
                });

                logger.trace("Reading mean expression");

                groupMeanExpression.forEach((group, meanFile) -> {
                    Map<String, Double> ensgMeanExpression;
                    try {
                        ensgMeanExpression = readLines(meanFile).stream().skip(1).map(line -> line.split("\t")).map(
                                split -> Pair.of(split[0], Double.parseDouble(split[1]))).collect(
                                toMap(Pair::getKey, Pair::getValue));
                    } catch (IOException e) {
                        throw new UncheckedIOException(e);
                    }
                    groups.forEach(tfGroup -> tfGroup.getTranscriptionFactors().forEach(
                            tf -> tf.setMeanExpression(group,
                                    tf.getEnsgs().stream().filter(ensgMeanExpression::containsKey).mapToDouble(
                                            ensgMeanExpression::get).average().orElse(-1))));
                });

                logger.trace("Reading tpm");

                groupTpm.forEach((group, tpmFile) -> {
                    Map<String, Double> ensgTpm;
                    try {
                        ensgTpm = readLines(tpmFile).stream().skip(1).map(line -> line.split("\t")).map(
                                split -> Pair.of(split[0], Double.parseDouble(split[1]))).collect(
                                toMap(Pair::getKey, Pair::getValue));
                    } catch (IOException e) {
                        throw new UncheckedIOException(e);
                    }
                    groups.forEach(tfGroup -> tfGroup.getTranscriptionFactors().forEach(tf -> tf.setTpm(group,
                            tf.getEnsgs().stream().filter(ensgTpm::containsKey).mapToDouble(
                                    ensgTpm::get).average().orElse(0.0))));
                });

                logger.trace("Reading ranks");

                hmRanks.forEach((hm, rankFile) -> {
                    Map<String, Integer> symbolRank;
                    try {
                        symbolRank = readLines(rankFile).stream().skip(1).map(line -> line.split("\t")).map(
                                split -> Pair.of(split[0], Integer.parseInt(split[1]))).collect(
                                toMap(Pair::getKey, Pair::getValue));
                    } catch (IOException e) {
                        throw new UncheckedIOException(e);
                    }
                    groups.forEach(group -> group.setRank(hm, symbolRank.getOrDefault(group.getSymbol(), -1)));
                });

                logger.trace("Writing output");

                JSONObject json = new JSONObject() {{
                    put("groups", groups.stream().map(TfGroup::toJSON).toList());
                }};

                File data = new File(outputFile, "data.json");

                try (FileWriter writer = new FileWriter(data)) {
                    json.write(writer, 4, 0);
                }

                return true;
            });
        }};
    }
}
