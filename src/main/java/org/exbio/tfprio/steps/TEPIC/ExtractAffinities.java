package org.exbio.tfprio.steps.TEPIC;

import org.apache.commons.lang3.tuple.Pair;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class ExtractAffinities extends ExecutableStep<Configs> {
    public final Map<String, Map<String, Collection<OutputFile>>> outputFiles = new HashMap<>();
    private final Map<InputFile, OutputFile> bridge = new HashMap<>();

    public ExtractAffinities(Configs configs, Map<String, Map<String, Collection<OutputFile>>> tepicFiles) {
        super(configs, false, tepicFiles.values().stream().flatMap(
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
            bridge.forEach((input, output) -> {
                add(() -> {
                    File[] allFiles = input.listFiles();

                    if (allFiles == null || allFiles.length == 0) {
                        return false;
                    }

                    List<String> attempts =
                            List.of("_Affinity_Gene_View_Filtered_TPM.txt", "_Affinity_Gene_View_Filtered.txt");

                    for (String attempt : attempts) {
                        List<File> matching =
                                Arrays.stream(allFiles).filter(file -> file.getName().endsWith(attempt)).toList();

                        if (matching.isEmpty()) {
                            continue;
                        }

                        if (matching.size() > 1) {
                            throw new IllegalStateException(
                                    "Found multiple files matching " + attempt + " in " + input);
                        }

                        File affinityFile = matching.get(0);

                        try (var reader = new BufferedReader(new FileReader(affinityFile));
                             var writer = new BufferedWriter(new FileWriter(output))) {
                            String header = reader.readLine().replaceAll("\\(.*?\\)", "");
                            String[] headerSplit = header.split("\t");

                            String first = headerSplit[0];
                            String[] tfs = Arrays.copyOfRange(headerSplit, 1, headerSplit.length);

                            List<String> tfOrder = Arrays.stream(tfs).distinct().sorted().toList();

                            writer.write(first + "\t" + String.join("\t", tfOrder));
                            writer.newLine();

                            reader.lines().map(line -> line.split("\t")).map(split -> {
                                String gene = split[0];
                                List<Double> affinities =
                                        Arrays.stream(split).skip(1).map(Double::parseDouble).toList();

                                Map<String, Double> tfAffinities = IntStream.range(0, affinities.size()).mapToObj(
                                        i -> Pair.of(tfs[i], affinities.get(i))).collect(
                                        Collectors.groupingBy(Pair::getKey,
                                                Collectors.averagingDouble(Pair::getValue)));

                                return Pair.of(gene, tfAffinities);
                            }).map(pair -> pair.getKey() + "\t" +
                                    tfOrder.stream().map(tf -> pair.getValue().getOrDefault(tf, 0.)).map(
                                            Object::toString).collect(Collectors.joining("\t"))).forEach(line -> {
                                try {
                                    writer.write(line);
                                    writer.newLine();
                                } catch (IOException e) {
                                    throw new UncheckedIOException(e);
                                }
                            });
                        }

                        return true;
                    }
                    logger.warn("Could not find affinity file for {}", input);
                    return false;
                });
            });
        }};
    }
}
