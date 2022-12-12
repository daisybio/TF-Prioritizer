package org.exbio.tfprio.steps.rnaSeq;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import static org.exbio.pipejar.util.FileManagement.readLines;
import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class DeSeq2 extends ExecutableStep {
    public final Collection<OutputFile> outputFiles = new HashSet<>();
    private final String scriptPath;
    private final InputFile batchFile;
    private final Map<InputFile, OutputFile> bridge;
    private final Map<InputFile, OutputFile> filteredBatchFiles = new HashMap<>();
    private final RequiredConfig<Integer> countFilter = new RequiredConfig<>(Configs.deSeq2.countFilter);

    public DeSeq2(Collection<OutputFile> dependencies, OutputFile batchFile) {
        super(false, dependencies, batchFile);

        bridge = new HashMap<>() {{
            dependencies.forEach(input -> {
                InputFile inputFile = addInput(input);
                OutputFile outputFile = addOutput(inputFile.getName());
                put(inputFile, outputFile);
                outputFiles.add(outputFile);
                filteredBatchFiles.put(inputFile, addOutput(
                        inputFile.getName().substring(0, inputFile.getName().lastIndexOf(".")) + "_meta.tsv"));
            });
        }};

        this.batchFile = addInput(batchFile);

        this.scriptPath = Objects.requireNonNull(getClass().getResource("deseq2.R")).getPath();
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            try {
                List<String> batchLines = readLines(batchFile);

                Map<String, String> sampleGroup = batchLines.stream().skip(1).map(line -> line.split("\t")).collect(
                        Collectors.toMap(line -> line[0], line -> line[1]));
                Map<String, String> sampleBatch = batchLines.stream().skip(1).map(line -> line.split("\t")).collect(
                        Collectors.toMap(line -> line[0], line -> line[2]));

                bridge.forEach((input, output) -> add(() -> {
                    OutputFile filteredBatchFile = filteredBatchFiles.get(input);

                    final List<String> order;
                    try (BufferedReader reader = new BufferedReader(new FileReader(input))) {
                        order = Arrays.stream(reader.readLine().split("\t")).skip(1).toList();
                    }

                    try (BufferedWriter writer = new BufferedWriter(new FileWriter(filteredBatchFile))) {
                        writer.write(batchLines.get(0));
                        writer.newLine();

                        order.forEach(sample -> {
                            try {
                                writer.write(
                                        String.join("\t", sample, sampleGroup.get(sample), sampleBatch.get(sample)));
                                writer.newLine();
                            } catch (IOException e) {
                                throw new RuntimeException(e);
                            }
                        });
                    }

                    executeAndWait("Rscript " + scriptPath + " --metadata " + filteredBatchFile.getAbsolutePath() +
                            " --input " + input.getAbsolutePath() + " --output " + output.getAbsolutePath() +
                            " --minCount " + countFilter.get(), true);
                    return true;
                }));
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }};
    }
}
