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
import java.util.stream.Stream;

import static org.exbio.pipejar.util.FileManagement.readLines;

public class Preprocessing extends ExecutableStep<Configs> {
    public final OutputFile outputFile = addOutput("preprocessed.tsv");

    private final Map<String, Pair<InputFile, InputFile>> inputFiles = new HashMap<>();

    public Preprocessing(Configs configs, Map<String, Pair<OutputFile, OutputFile>> hmThresholdGrouped) {
        super(configs, false, hmThresholdGrouped.values().stream().flatMap(
                stringOutputFileMap -> Stream.of(stringOutputFileMap.getLeft(),
                        stringOutputFileMap.getRight())).collect(Collectors.toSet()));

        hmThresholdGrouped.forEach((hm, outputFiles) -> {
            OutputFile inHm = new OutputFile(inputDirectory, hm);

            InputFile same = addInput(inHm, outputFiles.getLeft());
            InputFile different = addInput(inHm, outputFiles.getRight());
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
