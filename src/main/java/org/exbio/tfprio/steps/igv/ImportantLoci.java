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
import java.util.function.Function;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static org.exbio.pipejar.util.FileManagement.readLines;

public class ImportantLoci extends ExecutableStep<Configs> {
    public final Map<String, OutputFile> outputFiles = new HashMap<>();
    private final RequiredConfig<List<String>> importantLoci = new RequiredConfig<>(configs.igv.importantLoci);
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


    public ImportantLoci(Configs configs, OutputFile geneLocations, OutputFile ensgSymbol, OutputFile chipAtlas) {
        super(configs, false, geneLocations, ensgSymbol);

        this.ensgSymbolFile = addInput(ensgSymbol);
        this.geneInfo = addInput(geneLocations);
        this.chipAtlas = chipAtlas != null ? addInput(chipAtlas) : null;

        Set<String> groups = new HashSet<>();

        if (signalFiles.isSet()) {
            OutputFile signalDirectory = new OutputFile(inputDirectory, "signal");
            Arrays.stream(Objects.requireNonNull(signalFiles.get().listFiles(File::isDirectory))).forEach(groupDir -> {
                String group = groupDir.getName();
                OutputFile groupSignalDirectory = new OutputFile(signalDirectory, group);

                groups.add(group);

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

                        groups.add(group);

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

        if (groups.isEmpty()) {
            groups.add("all");
        }

        groups.forEach(group -> outputFiles.put(group, addOutput(group)));
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

            Map<String, Map<String, String>> patternDescriptionLocation =
                    importantLoci.get().stream().map(locusPattern -> {
                        Map<String, Set<String>> symbolEnsgs = ensgSymbol.entrySet().stream().filter(
                                entry -> entry.getValue().matches(locusPattern)).collect(
                                Collectors.groupingBy(Map.Entry::getValue,
                                        Collectors.mapping(Map.Entry::getKey, Collectors.toSet())));

                        Map<String, String> descriptionLocation = symbolEnsgs.entrySet().stream().flatMap(entry -> {
                            String symbol = entry.getKey();
                            Set<String> ensgs = entry.getValue();

                            return ensgs.stream().flatMap(ensg -> {
                                String description = symbol + (ensgs.size() > 1 ? " (" + ensg + ")" : "");

                                List<String> locations = geneLocations.get(ensg);

                                if (locations == null) {
                                    logger.warn("No location found for " + ensg + " (" + symbol + ")");
                                    return Stream.empty();
                                }

                                return IntStream.range(0, locations.size()).mapToObj(i -> Pair.of(description +
                                                (locations.size() > 1 ? " (" + (i + 1) + "/" + locations.size() + ")" : ""),
                                        locations.get(i)));
                            });
                        }).collect(Collectors.toMap(Pair::getKey, Pair::getValue));

                        return Pair.of(locusPattern, descriptionLocation);
                    }).collect(Collectors.toMap(Pair::getKey, Pair::getValue));

            outputFiles.forEach((group, outputFile) -> add(() -> {
                Collection<InputFile> experimentalFiles = groupExperimentalFiles.getOrDefault(group, new HashSet<>());
                Collection<InputFile> signalFiles = groupHmSignalFile.getOrDefault(group, new HashMap<>()).values();
                Function<File, String> removeExtension =
                        file -> file == null ? "" : file.getName().substring(0, file.getName().lastIndexOf('.'));
                Collection<File> allFiles = new HashSet<>() {{
                    addAll(experimentalFiles);
                    addAll(signalFiles);
                }};

                Map<File, String> descriptions = new HashMap<>(chipAtlasDescriptions) {{
                    experimentalFiles.forEach(file -> put(file, removeExtension.apply(file) + " (experimental)"));
                    signalFiles.forEach(file -> put(file, removeExtension.apply(file) + " (signal)"));
                }};

                IGV_Headless igv =
                        new IGV_Headless(logger, genome.get(), igvLib.get(), igvCacheDirectory.get(), outputFile);

                igv.createSession(allFiles.stream().toList(), chipAtlasFiles, descriptions);

                patternDescriptionLocation.forEach((pattern, descriptionLocation) -> {
                    File patternDirectory = new File(outputFile, pattern);
                    patternDirectory.mkdir();

                    igv.addCommand("snapshotDirectory " + patternDirectory.getAbsolutePath());

                    descriptionLocation.forEach((description, location) -> {
                        igv.addCommand("goto " + location);
                        igv.addCommand("snapshot " + description + ".png");
                    });
                });

                igv.run();

                return true;
            }));
        }};
    }
}
