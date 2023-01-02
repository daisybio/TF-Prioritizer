package org.exbio.tfprio.steps.distributionAnalysis;

import org.apache.commons.lang3.tuple.Pair;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.exbio.pipejar.util.FileManagement.readLines;

public class Preprocessing extends ExecutableStep {
    public final OutputFile outputFile = addOutput("preprocessed.tsv");

    private final Map<String, Pair<InputFile, InputFile>> inputFiles = new HashMap<>();

    public Preprocessing(Map<String, Map<Double, Pair<OutputFile, OutputFile>>> hmThresholdGrouped) {
        super(false, hmThresholdGrouped.values().stream().map(
                map -> map.get(map.keySet().stream().min(Double::compareTo).orElseThrow())).flatMap(
                maps -> Stream.of(maps.getLeft(), maps.getRight())).toList());

        hmThresholdGrouped.forEach((hm, thresholdMap) -> {
            OutputFile inHm = new OutputFile(inputDirectory, hm);
            var pair = thresholdMap.get(thresholdMap.keySet().stream().min(Double::compareTo).orElseThrow());

            InputFile same = addInput(inHm, pair.getLeft());
            InputFile different = addInput(inHm, pair.getRight());
            inputFiles.put(hm, Pair.of(same, different));
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                final Map<String, Map<String, HashSet<String>>> tfHmGrouped = new HashMap<>();
                inputFiles.forEach((hm, pair) -> {
                    InputFile same = pair.getLeft();
                    InputFile different = pair.getRight();

                    final Set<String> sameTfs;
                    final Set<String> differentTfs;
                    try {
                        sameTfs = getContainedTfs(same);
                        differentTfs = getContainedTfs(different);
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }

                    sameTfs.forEach(tf -> tfHmGrouped.computeIfAbsent(tf, x -> new HashMap<>()).computeIfAbsent(hm,
                            x -> new HashSet<>()).add("same"));
                    differentTfs.forEach(tf -> tfHmGrouped.computeIfAbsent(tf, x -> new HashMap<>()).computeIfAbsent(hm,
                            x -> new HashSet<>()).add("different"));
                });

                try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
                    writer.write("TF\tHM\tStages");
                    writer.newLine();

                    tfHmGrouped.forEach((tf, hmMap) -> {
                        hmMap.forEach((hm, stages) -> {
                            try {
                                writer.write(tf + "\t" + hm + "\t" + String.join(",", stages));
                                writer.newLine();
                            } catch (IOException e) {
                                throw new RuntimeException(e);
                            }
                        });
                    });
                }

                return true;
            });
        }};
    }

    private Set<String> getContainedTfs(InputFile input) throws IOException {
        return readLines(input).stream().skip(1).map(line -> line.split("\t")).map(line -> line[0]).collect(
                Collectors.toSet());
    }
}
