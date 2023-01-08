package org.exbio.tfprio.steps.rnaSeq;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.OptionalConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;

public class CreateBatchFile extends ExecutableStep {
    public final OutputFile outputFile = addOutput("batches.tsv");
    private final OptionalConfig<Map<String, Integer>> batches = new OptionalConfig<>(Configs.deSeq2.batches, false);
    private final Map<String, InputFile> files = new HashMap<>();

    public CreateBatchFile(Map<String, OutputFile> dependencies) {
        super(false, dependencies.values());

        dependencies.forEach((group, file) -> files.put(group, addInput(file)));
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                Map<String, Set<String>> samples = new HashMap<>() {{
                    files.forEach((group, inputFile) -> {
                        Set<String> currentSamples = new HashSet<>();
                        String cleanGroup = group.replace("-", "_");
                        put(cleanGroup, currentSamples);
                        try (BufferedReader reader = new BufferedReader(new FileReader(inputFile))) {
                            Arrays.stream(reader.readLine().split("\t")).filter(
                                    s -> !s.equals("gene_id") && !s.equals("Mean")).distinct().sorted().forEachOrdered(
                                    currentSamples::add);
                        } catch (IOException e) {
                            throw new RuntimeException(e);
                        }
                    });
                }};

                try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
                    writer.write("sample\tgroup\tbatch");
                    writer.newLine();
                    samples.forEach((group, sampleSet) -> sampleSet.forEach(sample -> {
                        if (batches.isSet() && !batches.get().containsKey(sample)) {
                            logger.warn("Sample {} is not in the batch file, using default batch 0", sample);
                        }

                        try {
                            writer.write(sample);
                            writer.write("\t");
                            writer.write(group);
                            writer.write("\t");
                            writer.write(batches.isSet() ? batches.get().getOrDefault(sample, 0).toString() : "1");
                            writer.newLine();
                        } catch (IOException e) {
                            throw new RuntimeException(e);
                        }
                    }));
                }

                return true;
            });
        }};
    }
}
