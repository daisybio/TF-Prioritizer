package org.exbio.tfprio.steps.distributionAnalysis;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.pipejar.util.FileFilters.Filters;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;

public class CreateBackground extends ExecutableStep {
    public final Map<String, OutputFile> outputFiles = new HashMap<>();
    private final Map<InputFile, OutputFile> bridge = new HashMap<>();

    public CreateBackground(Map<String, OutputFile> inputFiles) {
        super(false, inputFiles.values());

        inputFiles.forEach((hm, inDir) -> {
            OutputFile output = addOutput(hm + ".tsv");
            InputFile input = addInput(inDir);
            outputFiles.put(hm, output);
            bridge.put(input, output);
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((input, output) -> add(() -> {
                try (BufferedWriter writer = new BufferedWriter(new FileWriter(output))) {
                    writer.write("Gene\tTfTgScore\tHM\tPairing\tTF\tRegressionCoefficient");
                    writer.newLine();

                    Arrays.stream(Objects.requireNonNull(input.listFiles(Filters.fileFilter))).forEach(file -> {
                        try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
                            reader.lines().skip(1).forEachOrdered(line -> {
                                try {
                                    writer.write(line);
                                    writer.newLine();
                                } catch (IOException e) {
                                    throw new RuntimeException(e);
                                }
                            });
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
