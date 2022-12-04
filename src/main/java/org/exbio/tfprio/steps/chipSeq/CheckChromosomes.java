package org.exbio.tfprio.steps.chipSeq;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;

public class CheckChromosomes extends ExecutableStep {
    public final Map<String, Map<String, Collection<OutputFile>>> outputFiles = new HashMap<>();
    private final RequiredConfig<File> inputDirectory = new RequiredConfig<>(Configs.inputConfigs.chipSeq);
    private final Map<File, OutputFile> bridge = new HashMap<>();

    public CheckChromosomes() {
        super();

        // Register all input files and their corresponding output files
        InputFile input = addInput(inputDirectory);
        Arrays.stream(Objects.requireNonNull(input.listFiles())).forEach(d_group -> {
            OutputFile d_groupOut = addOutput(d_group.getName());
            d_groupOut.mkdir();
            outputFiles.put(d_group.getName(), new HashMap<>());
            Arrays.stream(Objects.requireNonNull(d_group.listFiles())).forEach(d_hm -> {
                OutputFile d_hmOut = addOutput(d_groupOut, d_hm.getName());
                d_hmOut.mkdir();
                outputFiles.get(d_group.getName()).put(d_hm.getName(), new HashSet<>());
                Arrays.stream(Objects.requireNonNull(d_hm.listFiles())).forEach(f_sample -> {
                    OutputFile f_sampleOut = addOutput(d_hmOut, f_sample.getName());
                    bridge.put(f_sample, f_sampleOut);
                    outputFiles.get(d_group.getName()).get(d_hm.getName()).add(f_sampleOut);
                });
            });
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((inputFile, outputFile) -> add(() -> {
                try (BufferedReader reader = new BufferedReader(new FileReader(inputFile));
                     BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
                    for (String line; (line = reader.readLine()) != null; ) {
                        // Remove initial "chr" from chromosome names
                        if (!line.startsWith("chr")) {
                            writer.write(line);
                        } else {
                            writer.write(line.substring("chr".length()));
                        }
                        writer.newLine();
                    }
                }
                return true;
            }));
        }};
    }
}
