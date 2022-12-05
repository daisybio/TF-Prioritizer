package org.exbio.tfprio.steps.tGene;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.*;
import java.util.Collection;
import java.util.HashSet;
import java.util.concurrent.Callable;

public class TGenePreprocess extends ExecutableStep {
    public final OutputFile outputFile;
    private final RequiredConfig<File> geneAnnotationFile =
            new RequiredConfig<>(Configs.inputConfigs.geneAnnotationFile);
    private final InputFile inputFile;

    public TGenePreprocess() {
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
                        if (line.startsWith("#") || !line.split("\t")[2].equals("transcript")) {
                            continue;
                        }

                        if (line.startsWith("chr")) {
                            line = line.substring("chr".length());
                        }

                        if (line.startsWith("M\t")) {
                            line = "MT\t" + line.substring("M\t".length());
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
