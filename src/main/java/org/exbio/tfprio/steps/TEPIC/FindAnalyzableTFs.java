package org.exbio.tfprio.steps.TEPIC;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;

public class FindAnalyzableTFs extends ExecutableStep {
    public final OutputFile outputFile = addOutput("analyzableTFs.tsv");

    public FindAnalyzableTFs(Map<String, Map<String, OutputFile>> affinityFiles) {
        super(true, affinityFiles.values().stream().flatMap(m -> m.values().stream()).toList());
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {

                List<String> tfs = getInputs().stream().flatMap(file -> {
                    try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
                        return Arrays.stream(reader.readLine().split("\t")).skip(1);
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                }).distinct().sorted().toList();

                try (var writer = new BufferedWriter(new FileWriter(outputFile))) {
                    tfs.forEach(tf -> {
                        try {
                            writer.write(tf);
                            writer.newLine();
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
