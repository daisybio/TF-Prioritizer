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

import static org.exbio.pipejar.util.FileManagement.extend;
import static org.exbio.pipejar.util.FileManagement.readLines;

public class DistributionTargetGenes extends ExecutableStep {
    public final Map<String, Map<String, OutputFile>> outputFiles = new HashMap<>();
    private final InputFile geneInfo;
    private final InputFile dcgFile;
    private final RequiredConfig<Integer> windowExtend =
            new RequiredConfig<>(org.exbio.tfprio.configs.Configs.igv.windowExtend);
    private final Map<InputFile, String> inputHm = new HashMap<>();
    private final Map<InputFile, String> inputGroup = new HashMap<>();
    private final OptionalConfig<File> experimentalDirectory =
            new OptionalConfig<>(Configs.igv.experimentalFiles, false);
    private final OptionalConfig<File> signalFiles = new OptionalConfig<>(Configs.igv.signalFiles, false);
    private final Map<String, Map<String, Collection<InputFile>>> groupHmTepicSampleDirectories = new HashMap<>();
    private final Map<String, Map<String, InputFile>> groupHmPredictedBindingSitesDirectories = new HashMap<>();
    private final Map<String, Set<InputFile>> groupExperimentalFiles = new HashMap<>();
    private final Map<String, Map<String, InputFile>> groupHmSignalFile = new HashMap<>();
    private final RequiredConfig<File> igvExecutable = new RequiredConfig<>(Configs.igv.igvExecutable);
    private final RequiredConfig<File> igvCacheDirectory = new RequiredConfig<>(Configs.igv.igvCacheDirectory);
    private final InputFile chipAtlas;
    private final InputFile ensgSymbolFile;
    private final RequiredConfig<String> genome = new RequiredConfig<>(Configs.igv.genome);
    private final Map<String, Map<String, InputFile>> groupHmTfTargetGenes = new HashMap<>();
    private final Map<Pair<String, String>, Collection<String>> groupPairingHms;

    public DistributionTargetGenes(OutputFile ensgSymbolFile, OutputFile geneInfo, OutputFile dcgFile,
                                   OutputFile chipAtlas, Map<String, Map<String, OutputFile>> groupHmTfTargetGenes,
                                   Map<String, Map<String, Collection<OutputFile>>> groupHmTepicSampleDirectories,
                                   Map<String, Map<String, OutputFile>> groupHmTepicBindingRegions) {
        super(false, groupHmTfTargetGenes.values().stream().flatMap(sub -> sub.values().stream()).toList(),
                groupHmTepicSampleDirectories.values().stream().flatMap(sub -> sub.values().stream()).flatMap(
                        Collection::stream).toList(),
                groupHmTepicBindingRegions.values().stream().flatMap(sub -> sub.values().stream()).toList(),
                new HashSet<>() {{
                    add(geneInfo);
                    add(dcgFile);
                    if (chipAtlas != null) {
                        add(chipAtlas);
                    }
                }});

        this.ensgSymbolFile = addInput(ensgSymbolFile);

        Set<String> groups = groupHmTfTargetGenes.keySet();
        Set<String> hms = groupHmTfTargetGenes.values().stream().flatMap(sub -> sub.keySet().stream()).collect(
                Collectors.toSet());

        groupPairingHms = groups.stream().flatMap(
                g1 -> groups.stream().filter(g2 -> g1.compareTo(g2) < 0).map(g2 -> Pair.of(g1, g2))).flatMap(
                pair -> hms.stream().filter(hm -> groupHmTfTargetGenes.get(pair.getLeft()).containsKey(hm) &&
                        groupHmTfTargetGenes.get(pair.getRight()).containsKey(hm)).map(
                        hm -> Pair.of(pair, hm))).collect(Collectors.groupingBy(Pair::getKey,
                Collectors.mapping(Pair::getValue, Collectors.toCollection(HashSet::new))));

        groupPairingHms.forEach((pair, availableHms) -> {
            String pairName = pair.getLeft() + "_" + pair.getRight();
            OutputFile outputPair = new OutputFile(outputDirectory, pairName);
            availableHms.forEach(hm -> {
                OutputFile outputHm = addOutput(outputPair, hm);
                outputFiles.computeIfAbsent(pairName, k -> new HashMap<>()).put(hm, outputHm);
            });
        });

        this.geneInfo = addInput(geneInfo);
        this.dcgFile = addInput(dcgFile);

        if (chipAtlas != null) {
            this.chipAtlas = addInput(chipAtlas);
        } else {
            this.chipAtlas = null;
        }

        OutputFile inTargetGenes = new OutputFile(inputDirectory, "targetGenes");
        OutputFile inTepic = new OutputFile(inputDirectory, "tepic");
        OutputFile inTfChipSeq = new OutputFile(inputDirectory, "tfChipSeq");
        OutputFile inTdf = new OutputFile(inputDirectory, "tdf");
        OutputFile inTfPredictedBindingSites = new OutputFile(inputDirectory, "tfPredictedBindingSites");

        groupHmTfTargetGenes.forEach((group, hmTargetGenes) -> {
            OutputFile inTargetGenesGroup = new OutputFile(inTargetGenes, group);
            OutputFile inTepicGroup = new OutputFile(inTepic, group);
            OutputFile inTfPredictedBindingSitesGroup = new OutputFile(inTfPredictedBindingSites, group);

            if (experimentalDirectory.isSet()) {
                File tfChipSeqGroupDirectory = new File(experimentalDirectory.get(), group);
                if (!tfChipSeqGroupDirectory.exists()) {
                    logger.warn("No TF ChIP-Seq data for group " + group);
                } else {
                    OutputFile inTfChipSeqGroup = new OutputFile(inTfChipSeq, group);

                    Set<InputFile> tfChipSeqFiles =
                            Arrays.stream(Objects.requireNonNull(tfChipSeqGroupDirectory.listFiles())).map(
                                    directory -> Arrays.stream(
                                            Objects.requireNonNull(directory.listFiles())).findFirst()).filter(
                                    Optional::isPresent).map(Optional::get).map(
                                    file -> new OutputFile(file.getAbsolutePath())).map(
                                    file -> addInput(inTfChipSeqGroup, file)).collect(Collectors.toSet());

                    groupExperimentalFiles.put(group, tfChipSeqFiles);
                }
            }

            hmTargetGenes.forEach((hm, targetGenes) -> {
                InputFile inTargetGenesHm = addInput(inTargetGenesGroup, targetGenes);
                OutputFile inTepicHm = new OutputFile(inTepicGroup, hm);
                InputFile inTfPredictedBindingSitesHm =
                        addInput(inTfPredictedBindingSitesGroup, groupHmTepicBindingRegions.get(group).get(hm));

                groupHmPredictedBindingSitesDirectories.computeIfAbsent(group, k -> new HashMap<>()).put(hm,
                        inTfPredictedBindingSitesHm);

                if (signalFiles.isSet()) {
                    File signalHmDirectory = extend(signalFiles.get(), group, hm);
                    if (signalHmDirectory.exists()) {
                        File tdfFile =
                                Arrays.stream(Objects.requireNonNull(signalHmDirectory.listFiles())).findFirst().orElse(
                                        null);
                        if (tdfFile != null) {
                            InputFile inTdfHm = addInput(inTdf, new OutputFile(tdfFile.getAbsolutePath()));
                            groupHmSignalFile.computeIfAbsent(group, k -> new HashMap<>()).put(hm, inTdfHm);
                        } else {
                            logger.warn("No signal file for group " + group + " and hm " + hm);
                        }
                    } else {
                        logger.warn("No signal directory for group " + group + " and hm " + hm);
                    }
                }

                groupHmTepicSampleDirectories.get(group).get(hm).forEach(sample -> {
                    InputFile inTepicSample = addInput(inTepicHm, sample);
                    this.groupHmTepicSampleDirectories.computeIfAbsent(group, g -> new HashMap<>()).computeIfAbsent(hm,
                            h -> new HashSet<>()).add(inTepicSample);
                });

                this.groupHmTfTargetGenes.computeIfAbsent(group, k -> new HashMap<>()).put(hm, inTargetGenesHm);
                inputHm.put(inTargetGenesHm, hm);
                inputGroup.put(inTargetGenesHm, group);
            });
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            final Map<String, List<String>> geneLocations;
            final Set<String> tfNames;
            final List<File> chipAtlasFiles;
            final Map<String, Map<String, Map<String, File>>> groupHmTfPredictedBindingSites;
            final Map<String, String> ensgSymbol;
            final Map<File, String> chipAtlasDescriptions;

            try {
                geneLocations = readLines(geneInfo).stream().map(line -> line.split("\t")).map(
                        split -> Pair.of(split[0],
                                "chr" + split[2] + ":" + (Integer.parseInt(split[3]) - windowExtend.get()) + "-" +
                                        (Integer.parseInt(split[4]) + windowExtend.get()))).distinct().collect(
                        Collectors.groupingBy(Pair::getKey, Collectors.mapping(Pair::getValue, Collectors.toList())));


                tfNames = readLines(dcgFile).stream().skip(1).map(line -> line.split("\t")[0]).collect(
                        Collectors.toSet());
                chipAtlasFiles = chipAtlas != null ? Arrays.stream(
                        Objects.requireNonNull(chipAtlas.listFiles(file -> file.getName().endsWith(".bed")))).toList() :
                        new ArrayList<>();
                groupHmTfPredictedBindingSites = groupHmPredictedBindingSitesDirectories.entrySet().stream().collect(
                        Collectors.toMap(Map.Entry::getKey,
                                groupEntry -> groupEntry.getValue().entrySet().stream().collect(
                                        Collectors.toMap(Map.Entry::getKey, hmEntry -> Arrays.stream(
                                                Objects.requireNonNull(hmEntry.getValue().listFiles(File::isFile))).map(
                                                tfFile -> {
                                                    String fileName = tfFile.getName();
                                                    String tfName = fileName.substring(0, fileName.lastIndexOf('.'));
                                                    return Pair.of(tfName, tfFile);
                                                }).collect(Collectors.toMap(Pair::getKey, Pair::getValue))))));
                ensgSymbol = readLines(ensgSymbolFile).stream().map(line -> line.split("\t")).collect(
                        Collectors.toMap(line -> line[0], line -> line.length > 1 ? line[1] : line[0]));

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

            groupPairingHms.forEach((pair, hms) -> {
                String group1 = pair.getLeft();
                String group2 = pair.getRight();

                Pair<Collection<InputFile>, Collection<InputFile>> experimentalFiles =
                        Pair.of(groupExperimentalFiles.get(group1), groupExperimentalFiles.get(group2));

                hms.forEach(hm -> {
                    Function<String, Collection<File>> groupToTepicBedFiles =
                            group -> groupHmTepicSampleDirectories.getOrDefault(group, new HashMap<>()).getOrDefault(hm,
                                    new HashSet<>()).stream().map(
                                    directory -> directory.listFiles(file -> file.getName().endsWith(".bed"))).filter(
                                    Objects::nonNull).map(
                                    files -> Arrays.stream(files).findFirst().orElse(null)).filter(
                                    Objects::nonNull).collect(Collectors.toSet());

                    Pair<Collection<File>, Collection<File>> tepicBedFiles =
                            Pair.of(groupToTepicBedFiles.apply(group1), groupToTepicBedFiles.apply(group2));

                    Pair<InputFile, InputFile> signalFiles =
                            Pair.of(groupHmSignalFile.getOrDefault(group1, new HashMap<>()).get(hm),
                                    groupHmSignalFile.getOrDefault(group2, new HashMap<>()).get(hm));

                    tfNames.forEach(tf -> add(() -> {
                        Function<String, List<String>> groupToTargetGenes = group -> {
                            try {
                                return readLines(
                                        new File(groupHmTfTargetGenes.get(group).get(hm), tf + ".tsv")).stream().skip(
                                        1).map(line -> line.split("\t")[0]).collect(Collectors.toList());
                            } catch (IOException e) {
                                throw new RuntimeException(e);
                            }
                        };

                        Pair<List<String>, List<String>> targetGenes =
                                Pair.of(groupToTargetGenes.apply(group1), groupToTargetGenes.apply(group2));

                        File tfOutputDir = new File(outputFiles.get(group1 + "_" + group2).get(hm), tf);

                        String query = tf.contains("::") ? tf + "_merged" : tf;

                        Pair<File, File> predictedBindingSites =
                                Pair.of(groupHmTfPredictedBindingSites.get(group1).get(hm).get(query),
                                        groupHmTfPredictedBindingSites.get(group2).get(hm).get(query));

                        Function<File, String> removeExtension = file -> file == null ? "" :
                                file.getName().substring(0, file.getName().lastIndexOf('.'));
                        Map<File, String> descriptions = new HashMap<>(chipAtlasDescriptions) {{
                            experimentalFiles.getLeft().forEach(
                                    file -> put(file, removeExtension.apply(file) + " experimental"));
                            experimentalFiles.getRight().forEach(
                                    file -> put(file, removeExtension.apply(file) + " experimental"));

                            tepicBedFiles.getLeft().forEach(
                                    file -> put(file, group1 + ", " + hm + " all predicted peaks"));
                            tepicBedFiles.getRight().forEach(
                                    file -> put(file, group2 + ", " + hm + " all predicted peaks"));

                            put(signalFiles.getLeft(), removeExtension.apply(signalFiles.getLeft()) + " signal");
                            put(signalFiles.getRight(), removeExtension.apply(signalFiles.getRight()) + " signal");

                            put(predictedBindingSites.getLeft(), group1 + ", " + hm + ", " + tf + " predicted");
                            put(predictedBindingSites.getRight(), group2 + ", " + hm + ", " + tf + " predicted");
                        }};

                        List<File> includedFiles = new ArrayList<File>() {{
                            addAll(experimentalFiles.getLeft());
                            addAll(tepicBedFiles.getLeft());
                            add(signalFiles.getLeft());
                            add(predictedBindingSites.getLeft());

                            addAll(experimentalFiles.getRight());
                            addAll(tepicBedFiles.getRight());
                            add(signalFiles.getRight());
                            add(predictedBindingSites.getRight());
                        }}.stream().filter(Objects::nonNull).filter(File::exists).toList();

                        IGV_Headless igv =
                                new IGV_Headless(logger, genome.get(), igvExecutable.get(), igvCacheDirectory.get(),
                                        tfOutputDir);
                        igv.createSession(includedFiles, chipAtlasFiles, descriptions);

                        igv.addCommand("snapshotDirectory " + tfOutputDir.getAbsolutePath());

                        Stream.of(targetGenes.getLeft(), targetGenes.getRight()).forEach(
                                targetGenesList -> IntStream.range(0, targetGenesList.size()).forEach(i -> {
                                    String targetGene = targetGenesList.get(i);

                                    if (!geneLocations.containsKey(targetGene)) {
                                        return;
                                    }

                                    List<String> locations = geneLocations.get(targetGene);

                                    IntStream.range(0, locations.size()).forEach(j -> {
                                        igv.addCommand("goto " + locations.get(j));
                                        igv.addCommand("snapshot " + (i + 1) + "_" + ensgSymbol.get(targetGene) +
                                                (locations.size() > 1 ? "[" + j + "]" : "") + ".png");
                                    });
                                }));

                        igv.run(tfOutputDir);

                        return true;
                    }));
                });
            });
        }};
    }
}
