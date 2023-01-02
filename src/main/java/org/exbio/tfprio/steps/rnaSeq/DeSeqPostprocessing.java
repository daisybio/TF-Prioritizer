package org.exbio.tfprio.steps.rnaSeq;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.Callable;

public class DeSeqPostprocessing extends ExecutableStep {
    public final Map<String, OutputFile> outputFiles = new HashMap<>();
    private final Map<InputFile, OutputFile> bridge;

    public DeSeqPostprocessing(Map<String, OutputFile> dependencies) {
        super(false, dependencies.values());

        bridge = new HashMap<>() {{
            dependencies.forEach((pairing, dependency) -> {
                InputFile inputFile = addInput(dependency);
                OutputFile outputFile = addOutput(inputFile.getName());
                put(inputFile, outputFile);
                outputFiles.put(pairing, outputFile);
            });
        }};
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((inputFile, outputFile) -> add(() -> {
                try (BufferedReader reader = new BufferedReader(new FileReader(inputFile));
                     BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
                    reader.readLine(); // skip header
                    writer.write("geneID\tlog2fc");
                    writer.newLine();

                    for (String line = reader.readLine(); line != null; line = reader.readLine()) {
                        String[] fields = line.split("\t");
                        writer.write(fields[0].replace("\"", "") + "\t" + fields[2]);
                        writer.newLine();
                    }
                }

                return true;
            }));
        }};
    }
}
