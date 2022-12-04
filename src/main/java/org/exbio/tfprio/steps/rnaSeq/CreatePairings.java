package org.exbio.tfprio.steps.rnaSeq;

import org.apache.commons.lang3.tuple.Pair;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.IntStream;

import static org.exbio.pipejar.util.FileManagement.readLines;

public class CreatePairings extends ExecutableStep {
    public final Map<String, OutputFile> outputFiles = new HashMap<>();
    private final Map<OutputFile, Pair<InputFile, InputFile>> bridge = new HashMap<>();

    public CreatePairings(Map<String, OutputFile> groupFiles) {
        super(false, groupFiles.values());

        Map<String, InputFile> inputFiles = groupFiles.entrySet().stream().collect(HashMap::new,
                (m, e) -> m.put(e.getKey(), addInput(e.getValue())), HashMap::putAll);

        inputFiles.forEach((group1, file1) -> inputFiles.entrySet().stream().filter(
                entry -> entry.getKey().compareTo(group1) > 0).forEach(entry -> {
            String group2 = entry.getKey();
            InputFile file2 = entry.getValue();

            String name = group1 + "_" + group2;
            OutputFile outputFile = addOutput(name + ".tsv");

            outputFiles.put(name, outputFile);

            bridge.put(outputFile, Pair.of(file1, file2));
        }));
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((outputFile, pair) -> add(() -> {
                List<String> leftLines = readLines(pair.getLeft());
                List<String> rightLines = readLines(pair.getRight());

                if (leftLines.size() != rightLines.size()) {
                    throw new IllegalArgumentException("The number of lines in the two files must be equal");
                }

                try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {

                    IntStream.range(0, leftLines.size()).mapToObj(i -> {
                        String right = rightLines.get(i);
                        String left = leftLines.get(i);
                        return left + "\t" + right.substring(right.indexOf('\t') + 1);
                    }).forEach(line -> {
                        try {
                            writer.write(line);
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
