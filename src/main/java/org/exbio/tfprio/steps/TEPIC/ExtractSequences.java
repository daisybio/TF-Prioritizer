package org.exbio.tfprio.steps.TEPIC;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;

public class ExtractSequences extends ExecutableStep {
    public final Map<String, Map<String, Collection<OutputFile>>> outputFiles = new HashMap<>();
    private final Map<InputFile, OutputFile> bridge = new HashMap<>();

    public ExtractSequences(Map<String, Map<String, Collection<OutputFile>>> tepicFiles) {
        super(false, tepicFiles.values().stream().flatMap(
                stringCollectionMap -> stringCollectionMap.values().stream().flatMap(Collection::stream)).toList());

        tepicFiles.forEach((group, hmMap) -> {
            OutputFile inputGroup = new OutputFile(inputDirectory, group);
            OutputFile outputGroup = new OutputFile(outputDirectory, group);

            hmMap.forEach((hm, sampleDirs) -> {
                OutputFile inputHm = new OutputFile(inputGroup, hm);
                OutputFile outputHm = new OutputFile(outputGroup, hm);

                sampleDirs.forEach(sampleDir -> {
                    InputFile inputDir = addInput(inputHm, sampleDir);
                    OutputFile outputFile = addOutput(outputHm, sampleDir.getName() + ".tsv");

                    bridge.put(inputDir, outputFile);
                    outputFiles.computeIfAbsent(group, s -> new HashMap<>()).computeIfAbsent(hm,
                            s -> new HashSet<>()).add(outputFile);
                });
            });
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((input, output) -> add(() -> {
                File[] allFiles = input.listFiles(file -> file.getName().equals("trap_sequences.tsv"));

                if (allFiles == null || allFiles.length == 0) {
                    logger.warn("No trap_sequences.tsv file found in " + input);
                    return false;
                }

                if (allFiles.length > 1) {
                    logger.warn("Multiple trap_sequences.tsv files found in " + input);
                    return false;
                }

                File sequenceFile = allFiles[0];

                try (var reader = new BufferedReader(new FileReader(sequenceFile));
                     var writer = new BufferedWriter(new FileWriter(output))) {
                    reader.lines().filter(line -> !line.startsWith("#")).map(line -> line.split("\t")).forEachOrdered(
                            split -> {
                                String tfWithVersion = split[0];
                                String tf = tfWithVersion.replaceAll("\\(.*?\\)", "");

                                try {
                                    writer.write(
                                            tf + "\t" + String.join("\t", Arrays.copyOfRange(split, 1, split.length)));
                                    writer.newLine();
                                } catch (IOException e) {
                                    throw new UncheckedIOException(e);
                                }
                            });
                }
                return true;
            }));
        }};
    }
}
