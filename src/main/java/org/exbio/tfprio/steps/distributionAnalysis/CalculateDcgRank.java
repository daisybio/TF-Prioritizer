package org.exbio.tfprio.steps.distributionAnalysis;

import org.apache.commons.lang3.tuple.Pair;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import static org.exbio.pipejar.util.FileManagement.readLines;

public class CalculateDcgRank extends ExecutableStep<Configs> {
    // ToDo: Use output of calculateDcgPerHm
    public final OutputFile outputFile = addOutput("dcg.tsv");
    private final Map<String, InputFile> hmRankFiles = new HashMap<>();

    public CalculateDcgRank(Configs configs, Map<String, OutputFile> hmRanks) {
        super(configs, false, hmRanks.values());

        hmRanks.forEach((hm, hmRankFile) -> hmRankFiles.put(hm, addInput(hmRankFile)));
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                Map<String, Map<String, Integer>> hmTfRank = hmRankFiles.entrySet().stream().map(hmRankFileEntry -> {
                    try {
                        return Pair.of(hmRankFileEntry.getKey(),
                                readLines(hmRankFileEntry.getValue()).stream().skip(1).map(line -> {
                                    String[] split = line.split("\t");
                                    return Pair.of(split[0], Integer.parseInt(split[1]));
                                }).collect(Collectors.toMap(Pair::getKey, Pair::getValue)));
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                }).collect(Collectors.toMap(Pair::getKey, Pair::getValue));

                Set<String> allTfs = hmTfRank.values().stream().flatMap(tfRank -> tfRank.keySet().stream()).collect(
                        Collectors.toSet());

                Map<String, Double> tfScore = allTfs.stream().map(tf -> {
                    double score = hmTfRank.values().stream().filter(tfRank -> tfRank.containsKey(tf)).mapToDouble(
                            tfRank -> (tfRank.size() - tfRank.get(tf)) / (double) tfRank.size()).sum();
                    return Pair.of(tf, score);
                }).collect(Collectors.toMap(Pair::getKey, Pair::getValue));

                try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
                    writer.write("TF\tSCORE");
                    writer.newLine();

                    tfScore.entrySet().stream().sorted(
                            (a, b) -> Double.compare(b.getValue(), a.getValue())).forEachOrdered(tfScoreEntry -> {
                        try {
                            writer.write(tfScoreEntry.getKey() + "\t" + tfScoreEntry.getValue());
                            writer.newLine();
                        } catch (IOException e) {
                            throw new RuntimeException(e);
                        }
                    });
                }

                return true;
            });
        }};
    }
}
