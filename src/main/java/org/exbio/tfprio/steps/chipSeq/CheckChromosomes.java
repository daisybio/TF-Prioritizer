package org.exbio.tfprio.steps.chipSeq;

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
import java.util.stream.Collectors;

public class CheckChromosomes extends ExecutableStep {
    public final Map<String, Map<String, Collection<OutputFile>>> outputFiles = new HashMap<>();
    private final Map<InputFile, OutputFile> bridge = new HashMap<>();

    public CheckChromosomes(Map<String, Map<String, Collection<OutputFile>>> peakFiles) {
        super(false, peakFiles.values().stream().flatMap(
                stringCollectionMap -> stringCollectionMap.values().stream()).flatMap(Collection::stream).collect(
                Collectors.toSet()));

        // manage output files
        peakFiles.forEach((group, hmMap) -> {
            outputFiles.put(group, new HashMap<>());
            OutputFile d_groupOut = addOutput(group);
            hmMap.forEach((hm, sampleFiles) -> {
                outputFiles.get(group).put(hm, new HashSet<>());
                OutputFile d_hmOut = addOutput(d_groupOut, hm);

                sampleFiles.forEach(sampleFile -> {
                    InputFile inputFile = addInput(sampleFile);

                    OutputFile outputFile = addOutput(d_hmOut, inputFile.getName());
                    outputFiles.get(group).get(hm).add(outputFile);
                    bridge.put(inputFile, outputFile);
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