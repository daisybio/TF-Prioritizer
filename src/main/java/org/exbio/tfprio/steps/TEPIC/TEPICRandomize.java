package org.exbio.tfprio.steps.TEPIC;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.OptionalConfig;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.exbio.pipejar.util.FileManagement.copyFile;

public class TEPICRandomize extends ExecutableStep<Configs> {
    public final Map<String, Map<String, Collection<OutputFile>>> outputFiles = new HashMap<>();
    private final Map<InputFile, OutputFile> bridge = new HashMap<>();
    private final RequiredConfig<String> sequenceFileName = new RequiredConfig<>(configs.tepic.sequenceFileName);
    private final OptionalConfig<Double> tpmFilter = new OptionalConfig<>(configs.deSeq2.tpmFilter, false);

    // TODO put this behind after extraction of affinity files

    public TEPICRandomize(Configs configs, Map<String, Map<String, Collection<OutputFile>>> tepicResults) {
        super(configs, false,
                tepicResults.values().stream().flatMap(x -> x.values().stream()).flatMap(Collection::stream).collect(
                        Collectors.toList()));
        tepicResults.forEach((group, hmMap) -> {
            outputFiles.put(group, new HashMap<>());
            OutputFile outputGroup = addOutput(group);
            hmMap.forEach((hm, sampleDirectories) -> {
                outputFiles.get(group).put(hm, new HashSet<>(sampleDirectories.size()));
                OutputFile outputHm = addOutput(outputGroup, hm);
                sampleDirectories.forEach(sampleDirectory -> {
                    InputFile inputDirectory = addInput(sampleDirectory);
                    OutputFile outputDirectory = addOutput(outputHm, sampleDirectory.getName());
                    bridge.put(inputDirectory, outputDirectory);
                    outputFiles.get(group).get(hm).add(outputDirectory);
                });
            });
        });
    }


    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((inputDirectory, outputDirectory) -> add(() -> {
                String suffix = tpmFilter.isSet() ? "_Gene_View_Filtered_TPM.txt" : "_Gene_View_Filtered.txt";

                File[] consideredFiles =
                        inputDirectory.listFiles(file -> file.isFile() && file.getName().endsWith(suffix));

                if (consideredFiles == null || consideredFiles.length == 0) {
                    logger.warn(
                            "No affinity file found in " + inputDirectory.getAbsolutePath() + "for suffix " + suffix);
                    return false;
                } else if (consideredFiles.length > 1) {
                    logger.warn("Multiple affinity files found in " + inputDirectory.getAbsolutePath() + "for suffix " +
                            suffix + ". Using first file.");
                    return false;
                }

                File affinities = consideredFiles[0];
                File sequences = new File(inputDirectory, sequenceFileName.get());

                if (!sequences.exists()) {
                    logger.warn("No " + sequenceFileName.get() + " found in " + inputDirectory);
                    return false;
                }

                Arrays.stream(Objects.requireNonNull(inputDirectory.listFiles(
                        file -> file.isFile() && !file.equals(affinities) && !file.equals(sequences)))).forEach(
                        inputFile -> {
                            try {
                                copyFile(inputFile, new File(outputDirectory, inputFile.getName()));
                            } catch (IOException e) {
                                throw new RuntimeException(e);
                            }
                        });


                File shuffledAffinities = new File(outputDirectory, affinities.getName());

                // Shuffle all columns of the affinity file except the first one
                try (BufferedReader reader = new BufferedReader(new FileReader(affinities));
                     BufferedWriter writer = new BufferedWriter(new FileWriter(shuffledAffinities))) {
                    String[] header = reader.readLine().split("\t");
                    List<Integer> columnOrder = IntStream.range(1, header.length).boxed().collect(Collectors.toList());
                    Collections.shuffle(columnOrder);
                    columnOrder.add(0, 0);

                    writer.write(columnOrder.stream().map(i -> header[i]).collect(Collectors.joining("\t")));

                    String line;
                    while ((line = reader.readLine()) != null) {
                        String[] split = line.split("\t");
                        String outputLine = columnOrder.stream().map(i -> split[i]).collect(Collectors.joining("\t"));
                        writer.write(outputLine);
                        writer.newLine();
                    }
                }

                File shuffledSequences = new File(outputDirectory, sequences.getName());

                // Get all TFs
                Set<String> tfNames = new HashSet<>();
                try (BufferedReader reader = new BufferedReader(new FileReader(sequences))) {
                    String line;
                    while ((line = reader.readLine()) != null) {
                        if (line.startsWith("#")) {
                            continue;
                        }
                        tfNames.add(line.split("\t")[0]);
                    }
                }

                // Shuffle all TFs
                try (BufferedReader reader = new BufferedReader(new FileReader(sequences));
                     BufferedWriter writer = new BufferedWriter(new FileWriter(shuffledSequences))) {
                    String line;
                    while ((line = reader.readLine()) != null) {
                        if (line.startsWith("#") || line.startsWith("TF\t")) {
                            writer.write(line);
                            writer.newLine();
                            continue;
                        }

                        String[] split = line.split("\t");
                        String randomTf =
                                tfNames.stream().skip((int) (tfNames.size() * Math.random())).findFirst().get();
                        split[0] = randomTf;
                        writer.write(String.join("\t", split));
                        writer.newLine();
                    }
                }

                return true;
            }));
        }};
    }
}
