package org.exbio.tfprio.steps.distributionAnalysis;

import org.apache.commons.lang3.tuple.Pair;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import static org.exbio.pipejar.util.FileManagement.readLines;

public class CalculateDcgPerHm extends ExecutableStep<Configs> {
    public final Map<String, OutputFile> outputFiles = new HashMap<>();
    private final Map<InputFile, OutputFile> bridge = new HashMap<>();

    public CalculateDcgPerHm(Configs configs, Map<String, OutputFile> hmRanks) {
        super(configs, false, hmRanks.values());

        hmRanks.forEach((hm, hmRankFile) -> {
            OutputFile outputFile = addOutput(hm + ".tsv");
            InputFile inputFile = addInput(hmRankFile);
            bridge.put(inputFile, outputFile);
            outputFiles.put(hm, outputFile);
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((inputFile, outputFile) -> add(() -> {
                Map<String, Integer> tfRank = readLines(inputFile).stream().skip(1).map(line -> {
                    String[] split = line.split("\t");
                    return Pair.of(split[0], Integer.parseInt(split[1]));
                }).collect(Collectors.toMap(Pair::getKey, Pair::getValue));

                Map<String, Double> tfDcg = tfRank.entrySet().stream().map(tfRankEntry -> {
                    double dcg = (tfRank.size() - tfRankEntry.getValue()) / (double) tfRank.size();
                    return Pair.of(tfRankEntry.getKey(), dcg);
                }).collect(Collectors.toMap(Pair::getKey, Pair::getValue));

                try (var writer = new BufferedWriter(new FileWriter(outputFile))) {
                    tfDcg.forEach((tf, dcg) -> {
                        try {
                            writer.write(tf + "\t" + dcg);
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
