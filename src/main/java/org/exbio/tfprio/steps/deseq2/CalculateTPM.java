package org.exbio.tfprio.steps.deseq2;

import org.apache.commons.lang3.tuple.Pair;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.Callable;

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
                    Map<String, Double> meanCounts =
                            readLines(inputFile).stream().skip(1).collect(HashMap::new, (map, line) -> {
                                String[] split = line.split("\t");
                                map.put(split[0], Double.parseDouble(split[1]));
                            }, HashMap::putAll);


                    Map<String, Double> intermediate =
                            meanCounts.entrySet().stream().filter(entry -> lengths.containsKey(entry.getKey())).collect(
                                    HashMap::new, (map, entry) -> map.put(entry.getKey(),
                                            1E3 * entry.getValue() / lengths.get(entry.getKey())), HashMap::putAll);

                    double sumIntermediate = intermediate.values().stream().mapToDouble(Double::doubleValue).sum();

                    try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
                        intermediate.entrySet().stream().map(entry -> Pair.of(entry.getKey(),
                                (entry.getValue() / sumIntermediate) * 1E6)).sorted().forEachOrdered(pair -> {
                            try {
                                writer.write(pair.getLeft() + "\t" + pair.getRight());
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
