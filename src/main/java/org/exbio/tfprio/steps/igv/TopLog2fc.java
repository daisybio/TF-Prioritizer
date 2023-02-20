package org.exbio.tfprio.steps.igv;

import org.apache.commons.lang3.tuple.Pair;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.OptionalConfig;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;
import org.exbio.tfprio.lib.IGV_Headless;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static org.exbio.pipejar.util.FileManagement.readLines;

public class TopLog2fc extends ExecutableStep<Configs> {
    public final Map<String, Pair<OutputFile, OutputFile>> outputFiles = new HashMap<>();
    private final Map<String, InputFile> pairingDeseq2 = new HashMap<>();
    private final RequiredConfig<Integer> topLog2FoldChange = new RequiredConfig<>(configs.igv.topLog2FoldChange);
    private final InputFile ensgSymbolFile;
    private final RequiredConfig<File> igvLib = new RequiredConfig<>(configs.igv.igvLib);
    private final RequiredConfig<File> igvCacheDirectory = new RequiredConfig<>(configs.igv.igvCacheDirectory);
    private final InputFile geneInfo;
    private final RequiredConfig<Integer> windowExtend = new RequiredConfig<>(configs.igv.windowExtend);
    private final InputFile chipAtlas;
    private final RequiredConfig<String> genome = new RequiredConfig<>(configs.inputConfigs.genome);
    private final OptionalConfig<File> signalFiles = new OptionalConfig<>(configs.igv.signalFiles, false);
    private final OptionalConfig<File> experimentalDirectory =
            new OptionalConfig<>(configs.igv.experimentalFiles, false);
    private final Map<String, Map<String, InputFile>> groupHmSignalFile = new HashMap<>();
    private final Map<String, Set<InputFile>> groupExperimentalFiles = new HashMap<>();


    public TopLog2fc(Configs configs, Map<String, OutputFile> pairingDeseq2, OutputFile geneLocations,
                     OutputFile ensgSymbol, OutputFile chipAtlas) {
        super(configs, false, pairingDeseq2.values(), geneLocations, ensgSymbol, chipAtlas);

        this.ensgSymbolFile = addInput(ensgSymbol);
        this.geneInfo = addInput(geneLocations);
        this.chipAtlas = chipAtlas != null ? addInput(chipAtlas) : null;


        OutputFile deseq2 = new OutputFile(inputDirectory, "deseq2");
        pairingDeseq2.forEach((pairing, output) -> this.pairingDeseq2.put(pairing, addInput(deseq2, output)));

        if (signalFiles.isSet()) {
            OutputFile signalDirectory = new OutputFile(inputDirectory, "signal");
            Arrays.stream(Objects.requireNonNull(signalFiles.get().listFiles(File::isDirectory))).forEach(groupDir -> {
                String group = groupDir.getName();
                OutputFile groupSignalDirectory = new OutputFile(signalDirectory, group);

                File[] hmDirs = groupDir.listFiles();

                if (hmDirs == null) {
                    return;
                }

                Arrays.stream(hmDirs).forEach(hmDir -> {
                    String hm = hmDir.getName();
                    Arrays.stream(Objects.requireNonNull(hmDir.listFiles())).findFirst().ifPresent(
                            tdfFile -> groupHmSignalFile.computeIfAbsent(group, k -> new HashMap<>()).put(hm,
                                    addInput(groupSignalDirectory, new OutputFile(tdfFile.getAbsolutePath()))));
                });
            });
        }

        if (experimentalDirectory.isSet()) {
            OutputFile experimentalSignalDirectory = new OutputFile(inputDirectory, "experimental");
            Arrays.stream(Objects.requireNonNull(experimentalDirectory.get().listFiles(File::isDirectory))).forEach(
                    groupDir -> {
                        String group = groupDir.getName();
                        OutputFile groupSignalDirectory = new OutputFile(experimentalSignalDirectory, group);

                        File[] tfDirs = groupDir.listFiles();

                        if (tfDirs == null) {
                            return;
                        }

                        Arrays.stream(tfDirs).map(tfDirectory -> Arrays.stream(
                                Objects.requireNonNull(tfDirectory.listFiles())).findFirst()).filter(
                                Optional::isPresent).map(Optional::get).forEach(
                                experimentalFile -> groupExperimentalFiles.computeIfAbsent(group,
                                        k -> new HashSet<>()).add(addInput(groupSignalDirectory,
                                        new OutputFile(experimentalFile.getAbsolutePath()))));
                    });
        }

        this.pairingDeseq2.keySet().forEach(pairing -> {
            OutputFile pairingDir = new OutputFile(outputDirectory, pairing);

            this.outputFiles.put(pairing,
                    Pair.of(addOutput(pairingDir, "upregulated"), addOutput(pairingDir, "downregulated")));
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            final Map<String, List<String>> geneLocations;
            final Map<String, String> ensgSymbol;
            final List<File> chipAtlasFiles;
            final Map<File, String> chipAtlasDescriptions;


            try {
                geneLocations = readLines(geneInfo).stream().map(line -> line.split("\t")).map(
                        split -> Pair.of(split[0],
                                "chr" + split[2] + ":" + (Integer.parseInt(split[3]) - windowExtend.get()) + "-" +
                                        (Integer.parseInt(split[4]) + windowExtend.get()))).distinct().collect(
                        Collectors.groupingBy(Pair::getKey, Collectors.mapping(Pair::getValue, Collectors.toList())));
                ensgSymbol = readLines(ensgSymbolFile).stream().map(line -> line.split("\t")).collect(
                        Collectors.toMap(line -> line[0], line -> line.length > 1 ? line[1] : line[0]));
                chipAtlasFiles = chipAtlas != null ? Arrays.stream(
                        Objects.requireNonNull(chipAtlas.listFiles(file -> file.getName().endsWith(".bed")))).toList() :
                        new ArrayList<>();
                chipAtlasDescriptions = chipAtlasFiles.stream().collect(Collectors.toMap(file -> file, file -> {
                    try (var reader = new BufferedReader(new FileReader(file))) {
                        String header = reader.readLine();
                        String pattern = "name=\"(.*?)\"";
                        Pattern r = Pattern.compile(pattern);
                        Matcher m = r.matcher(header);

                        if (m.find()) {
                            return m.group(1);
                        } else {
                            return file.getName();
                        }
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                }));
            } catch (IOException e) {
                throw new RuntimeException(e);
            }

            pairingDeseq2.forEach((pairing, deseq2) -> add(() -> {
                List<String> sortedGenes = readLines(deseq2).stream().skip(1).map(line -> {
                    String[] split = line.split("\t");

                    return Pair.of(split[0].replace("\"", ""), Double.parseDouble(split[2]));
                }).sorted(Comparator.comparingDouble(Pair::getRight)).map(Pair::getKey).toList();

                List<String> upregulated = sortedGenes.subList(0, topLog2FoldChange.get()).stream().toList();
                List<String> downregulated = sortedGenes.subList(sortedGenes.size() - topLog2FoldChange.get(),
                        sortedGenes.size()).stream().toList();

                // Create a Map<String, String> from gene symbol to gene location based on the interesting gene IDs.
                // If multiple IDs are associated with the same gene symbol, add a counter to the gene symbol.

                Function<Collection<String>, Map<String, String>> getGeneSymbolLocations =
                        input -> input.stream().flatMap(geneID -> {
                            String symbol = ensgSymbol.getOrDefault(geneID, geneID);
                            return geneLocations.getOrDefault(geneID, new ArrayList<>()).stream().map(
                                    location -> Pair.of(symbol, location));
                        }).collect(Collectors.groupingBy(Pair::getKey,
                                Collectors.mapping(Pair::getValue, Collectors.toList()))).entrySet().stream().flatMap(
                                entry -> {
                                    if (entry.getValue().size() == 1) {
                                        return Stream.of(Pair.of(entry.getKey(), entry.getValue().get(0)));
                                    } else {
                                        return IntStream.range(0, entry.getValue().size()).mapToObj(
                                                i -> Pair.of(entry.getKey() + "_" + i, entry.getValue().get(i)));
                                    }
                                }).collect(Collectors.toMap(Pair::getKey, Pair::getValue));

                Map<String, String> upregulatedGeneSymbolLocations = getGeneSymbolLocations.apply(upregulated);
                Map<String, String> downregulatedGeneSymbolLocations = getGeneSymbolLocations.apply(downregulated);

                String group1 = pairing.split("_")[0];
                String group2 = pairing.split("_")[1];

                Collection<? extends File> group1SignalFiles = groupHmSignalFile.get(group1).values();
                Collection<? extends File> group2SignalFiles = groupHmSignalFile.get(group2).values();

                Collection<? extends File> group1ExperimentalFiles = groupExperimentalFiles.get(group1);
                Collection<? extends File> group2ExperimentalFiles = groupExperimentalFiles.get(group2);

                OutputFile upregulatedDir = outputFiles.get(pairing).getLeft();
                OutputFile downregulatedDir = outputFiles.get(pairing).getRight();

                List<File> allFiles = new ArrayList<>() {{
                    addAll(group1SignalFiles);
                    addAll(group1ExperimentalFiles);
                    addAll(group2SignalFiles);
                    addAll(group2ExperimentalFiles);
                }};
                Function<File, String> removeExtension =
                        file -> file == null ? "" : file.getName().substring(0, file.getName().lastIndexOf('.'));

                Map<File, String> descriptions = new HashMap<>(chipAtlasDescriptions) {{
                    group1ExperimentalFiles.forEach(file -> put(file, removeExtension.apply(file) + " (experimental)"));
                    group1SignalFiles.forEach(file -> put(file, removeExtension.apply(file) + " (signal)"));
                    group2ExperimentalFiles.forEach(file -> put(file, removeExtension.apply(file) + " (experimental)"));
                    group2SignalFiles.forEach(file -> put(file, removeExtension.apply(file) + " (signal)"));
                }};

                IGV_Headless igv = new IGV_Headless(logger, genome.get(), igvLib.get(), igvCacheDirectory.get(),
                        upregulatedDir.getParentFile());

                igv.createSession(allFiles, chipAtlasFiles, descriptions);

                Consumer<Map<String, String>> takeSnapshots =
                        geneSymbolLocations -> geneSymbolLocations.forEach((geneSymbol, location) -> {
                            igv.addCommand("goto " + location);
                            igv.addCommand("snapshot " + geneSymbol + ".png");
                        });

                igv.addCommand("snapshotDirectory " + upregulatedDir.getAbsolutePath());
                takeSnapshots.accept(upregulatedGeneSymbolLocations);

                igv.addCommand("snapshotDirectory " + downregulatedDir.getAbsolutePath());
                takeSnapshots.accept(downregulatedGeneSymbolLocations);

                igv.run();

                return true;
            }));
        }};
    }
}
