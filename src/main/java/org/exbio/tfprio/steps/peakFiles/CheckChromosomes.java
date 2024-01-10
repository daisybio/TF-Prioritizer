package org.exbio.tfprio.steps.peakFiles;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

public class CheckChromosomes extends ExecutableStep<Configs> {
    public final Map<String, Map<String, Collection<OutputFile>>> outputFiles = new HashMap<>();
    private final Map<InputFile, OutputFile> bridge = new HashMap<>();

    public CheckChromosomes(Configs configs, Map<String, Map<String, Collection<OutputFile>>> peakFiles) {
        super(configs, false, peakFiles.values().stream().flatMap(
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
                        String[] fields = line.split("\t");
                        String chromosome = fields[0];
                        String[] rest = Arrays.copyOfRange(fields, 1, fields.length);

                        // Remove "chr" prefix
                        if (chromosome.startsWith("chr")) {
                            chromosome = chromosome.substring("chr".length());
                        }

                        chromosome = chromosome.toUpperCase();

                        // Skip everything but numeric, X, Y, MT and M chromosomes
                        if (!chromosome.matches("^(\\d+|X|Y|MT|M)$")) {
                            continue;
                        }

                        String correctedLine = chromosome + "\t" + String.join("\t", rest);
                        writer.write(correctedLine);
                        
                        writer.newLine();
                    }
                }
                return true;
            }));
        }};
    }
}