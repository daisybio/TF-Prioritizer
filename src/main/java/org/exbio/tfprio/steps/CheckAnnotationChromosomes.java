package org.exbio.tfprio.steps;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.*;
import java.util.Collection;
import java.util.HashSet;
import java.util.concurrent.Callable;

public class CheckAnnotationChromosomes extends ExecutableStep {
    public final OutputFile outputFile;
    private final RequiredConfig<File> geneAnnotationFile =
            new RequiredConfig<>(Configs.inputConfigs.geneAnnotationFile);
    private final InputFile inputFile;

    public CheckAnnotationChromosomes() {
        super();

        inputFile = addInput(geneAnnotationFile);
        outputFile = addOutput("annotations.tsv");
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                try (BufferedReader reader = new BufferedReader(new FileReader(inputFile));
                     BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
                    for (String line; (line = reader.readLine()) != null; ) {
                        // Remove initial "chr" from chromosome names
                        if (line.startsWith("#")) {
                            continue;
                        }

                        if (line.startsWith("chr")) {
                            line = line.substring("chr".length());
                        }

                        writer.write(line);
                        writer.newLine();
                    }
                }
                return true;
            });
        }};
    }
}
