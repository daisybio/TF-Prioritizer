package org.exbio.tfprio.steps.rnaSeq;

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
import java.util.stream.IntStream;

import static org.exbio.pipejar.util.FileManagement.readLines;


public class CalculateTPM extends ExecutableStep<Configs> {
    public final Map<String, OutputFile> outputFiles = new HashMap<>();
    private final InputFile lengthsFile;
    private final Map<InputFile, OutputFile> bridge = new HashMap<>();

    public CalculateTPM(Configs configs, Map<String, OutputFile> countFiles, OutputFile lengthsFile) {
        super(configs, false, countFiles.values(), lengthsFile);
        this.lengthsFile = addInput(lengthsFile);

        countFiles.forEach((group, input) -> {
            InputFile inputFile = addInput(input);
            OutputFile outputFile = addOutput(inputFile.getName());
            outputFiles.put(group, outputFile);
            bridge.put(inputFile, outputFile);
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            Map<String, Integer> lengths;
            try {
                lengths = readLines(lengthsFile).stream().collect(HashMap::new, (map, line) -> {
                    String[] split = line.split("\t");
                    if (!split[1].equals("null")) {
                        map.put(split[0], Integer.parseInt(split[1]));
                    }
                }, HashMap::putAll);
            } catch (IOException e) {
                throw new RuntimeException(e);
            }


            bridge.forEach((inputFile, outputFile) -> add(() -> {
                List<String> inputLines = readLines(inputFile);
                Map<String, List<Double>> counts = inputLines.stream().skip(1).collect(HashMap::new, (map, line) -> {
                    String[] split = line.split("\t");
                    map.put(split[0], Arrays.stream(split).skip(1).map(Double::parseDouble).toList());
                }, HashMap::putAll);

                logger.trace("Found " + counts.size() + " lines.");

                Map<String, List<Double>> intermediate =
                        counts.entrySet().stream().filter(entry -> lengths.containsKey(entry.getKey())).collect(
                                Collectors.toMap(Map.Entry::getKey,
                                        entry -> entry.getValue().stream().mapToDouble(Double::doubleValue).map(
                                                value -> 1E3 * value / lengths.get(entry.getKey())).boxed().collect(
                                                Collectors.toList())));

                logger.trace("Found " + intermediate.size() + " lines with known lengths.");

                List<Double> sumIntermediate =
                        IntStream.range(0, intermediate.values().iterator().next().size()).mapToObj(
                                i -> intermediate.values().stream().mapToDouble(list -> list.get(i)).sum()).toList();

                try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
                    writer.write(inputLines.get(0));
                    writer.newLine();
                    intermediate.entrySet().stream().map(entry -> entry.getKey() + "\t" +
                            IntStream.range(0, entry.getValue().size()).mapToObj(
                                    i -> (entry.getValue().get(i) / sumIntermediate.get(i) * 1E6)).map(
                                    String::valueOf).collect(Collectors.joining("\t"))).sorted().forEachOrdered(
                            line -> {
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
