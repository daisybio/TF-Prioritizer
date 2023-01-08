package org.exbio.tfprio.steps.rnaSeq;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.*;
import java.util.Collection;
import java.util.HashSet;
import java.util.concurrent.Callable;

public class FilterEnsgs extends ExecutableStep {
    public final OutputFile outputFile = addOutput("ensgs.tsv");
    public final OutputFile cleanFile = addOutput("clean.tsv");
    private final RequiredConfig<File> geneIdsFile = new RequiredConfig<>(Configs.inputConfigs.geneIDs);
    private final InputFile inputFile = addInput(geneIdsFile);

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                try (var writer = new BufferedWriter(new FileWriter(outputFile));
                     var cleanWriter = new BufferedWriter(new FileWriter(cleanFile));
                     var reader = new BufferedReader(new FileReader(inputFile))) {
                    reader.lines().forEach(id -> {
                        try {
                            if (id.startsWith("ENSMUSG") || id.startsWith("ENSG")) {
                                String withoutVersion = id.split("\\.")[0];
                                writer.write(withoutVersion);
                                writer.newLine();

                                cleanWriter.write(withoutVersion);
                                cleanWriter.newLine();
                            } else {
                                writer.write("-");
                                writer.newLine();
                            }
                        } catch (IOException e) {
                            throw new RuntimeException(e);
                        }
                    });
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }

                return true;
            });
        }};
    }
}
