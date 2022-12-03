package org.exbio.tfprio.steps.deseq2;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.FileManagement.readLines;

public class FilterTPM extends ExecutableStep {
    private final HashMap<InputFile, OutputFile> bridge;

    private final RequiredConfig<Double> tpmFilter = new RequiredConfig<>(Configs.deSeq2.tpmFilter);

    public FilterTPM(Collection<OutputFile> dependencies) {
        super(false, dependencies);

        bridge = dependencies.stream().collect(HashMap::new,
                (map, dependency) -> map.put(addInput(dependency), addOutput(dependency.getName())), HashMap::putAll);
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((inputFile, outputFile) -> add(() -> {
                try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
                    readLines(inputFile).stream().map(line -> line.split("\t")).forEachOrdered(split -> {
                        try {
                            writer.write(String.join("\t", split[0],
                                    Double.parseDouble(split[1]) >= tpmFilter.get() ? split[1] : "0"));
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
