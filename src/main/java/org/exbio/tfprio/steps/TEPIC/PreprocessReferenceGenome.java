package org.exbio.tfprio.steps.TEPIC;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.*;
import java.util.Collection;
import java.util.HashSet;
import java.util.concurrent.Callable;

public class PreprocessReferenceGenome extends ExecutableStep<Configs> {
    public final OutputFile outputFile = addOutput("referenceGenome.fa");
    private final RequiredConfig<File> referenceGenomeConfig = new RequiredConfig<>(configs.tepic.inputReferenceGenome);
    private final InputFile referenceGenome = addInput(referenceGenomeConfig);

    public PreprocessReferenceGenome(Configs configs) {
        super(configs);
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                try (var reader = new BufferedReader(new FileReader(referenceGenome));
                     var writer = new BufferedWriter(new FileWriter(outputFile))) {
                    reader.lines().map(line -> line.startsWith(">") ? line.replace("chr", "") : line).forEach(line -> {
                        try {
                            writer.write(line);
                            writer.newLine();
                        } catch (IOException e) {
                            throw new RuntimeException(e);
                        }
                    });
                }
                return true;
            });
        }};
    }
}
