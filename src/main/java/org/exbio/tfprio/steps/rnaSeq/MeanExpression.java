package org.exbio.tfprio.steps.rnaSeq;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.*;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.Callable;

public class MeanExpression extends ExecutableStep<Configs> {
    public final Map<String, OutputFile> outputFiles = new HashMap<>();
    private final Map<InputFile, OutputFile> bridge = new HashMap<>();

    public MeanExpression(Configs configs, Map<String, OutputFile> inputFiles) {
        super(configs, false, inputFiles.values());

        inputFiles.forEach((group, input) -> {
            InputFile inputFile = addInput(input);
            OutputFile outputFile = addOutput(inputFile.getName());
            bridge.put(inputFile, outputFile);
            outputFiles.put(group, outputFile);
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((input, output) -> add(() -> {
                try (BufferedReader reader = new BufferedReader(new FileReader(input));
                     BufferedWriter writer = new BufferedWriter(new FileWriter(output))) {
                    reader.readLine(); // Skip header

                    writer.write("gene_id\tmean");
                    writer.newLine();

                    String line;
                    while ((line = reader.readLine()) != null) {
                        String[] split = line.split("\t");
                        double sum = 0;
                        for (int i = 1; i < split.length; i++) {
                            sum += Double.parseDouble(split[i]);
                        }
                        writer.write(split[0] + "\t" + sum / (split.length - 1));
                        writer.newLine();
                    }
                } catch (IOException e) {
                    throw new RuntimeException("Could not read/write file", e);
                }
                return true;
            }));
        }};
    }
}
