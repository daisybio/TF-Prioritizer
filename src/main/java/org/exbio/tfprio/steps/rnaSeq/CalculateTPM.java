package org.exbio.tfprio.steps.rnaSeq;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.exbio.pipejar.util.FileManagement.readLines;


public class CalculateTPM extends ExecutableStep {
    private final InputFile lengthsFile;
    private final Map<InputFile, OutputFile> bridge;

    public CalculateTPM(Collection<OutputFile> countFiles, OutputFile lengthsFile) {
        super(false, countFiles, lengthsFile);
        this.lengthsFile = addInput(lengthsFile);

        bridge = countFiles.stream().collect(HashMap::new,
                (map, dependency) -> map.put(addInput(dependency), addOutput(dependency.getName())), HashMap::putAll);
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            try {
                Map<String, Integer> lengths = readLines(lengthsFile).stream().collect(HashMap::new, (map, line) -> {
                    String[] split = line.split("\t");
                    if (!split[1].equals("null")) {
                        map.put(split[0], Integer.parseInt(split[1]));
                    }
                }, HashMap::putAll);

                bridge.forEach((inputFile, outputFile) -> add(() -> {
                    List<String> inputLines = readLines(inputFile);
                    Map<String, List<Double>> counts =
                            inputLines.stream().skip(1).collect(HashMap::new, (map, line) -> {
                                String[] split = line.split("\t");
                                map.put(split[0],
                                        Arrays.stream(split).skip(1).map(Double::parseDouble).collect(ArrayList::new,
                                                List::add, List::addAll));
                            }, HashMap::putAll);


                    Map<String, List<Double>> intermediate =
                            counts.entrySet().stream().filter(entry -> lengths.containsKey(entry.getKey())).collect(
                                    Collectors.toMap(Map.Entry::getKey,
                                            entry -> entry.getValue().stream().mapToDouble(Double::doubleValue).map(
                                                    value -> 1E3 * value / lengths.get(entry.getKey())).boxed().collect(
                                                    Collectors.toList())));

                    List<Double> sumIntermediate =
                            IntStream.range(0, intermediate.values().iterator().next().size()).mapToObj(
                                    i -> intermediate.values().stream().mapToDouble(
                                            list -> list.get(i)).sum()).toList();

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
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }};
    }

    @Override
    protected boolean mayBeSkipped() {
        return false;
    }
}
