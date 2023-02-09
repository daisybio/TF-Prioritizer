package org.exbio.tfprio.steps.chipAtlas;

import org.apache.commons.lang3.tuple.Pair;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;
import org.exbio.tfprio.lib.RegionWithPayload;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import static org.exbio.pipejar.util.FileManagement.makeSureFileExists;
import static org.exbio.pipejar.util.FileManagement.readLines;

public class CoOccurrenceBindingEnergies extends ExecutableStep<Configs> {
    public final Map<String, Map<String, OutputFile>> outputFiles = new HashMap<>();
    private final InputFile coOccurringRegions;
    private final Map<OutputFile, Collection<InputFile>> bridge = new HashMap<>();

    public CoOccurrenceBindingEnergies(Configs configs, Map<String, Map<String, Collection<OutputFile>>> affinities,
                                       OutputFile coOccurringRegions) {
        super(configs, false,
                affinities.values().stream().flatMap(m -> m.values().stream().flatMap(Collection::stream)).toList(),
                coOccurringRegions);

        this.coOccurringRegions = addInput(coOccurringRegions);

        affinities.forEach((group, hmAffinities) -> {
            OutputFile groupInput = new OutputFile(inputDirectory, group);
            OutputFile groupOutput = new OutputFile(outputDirectory, group);

            hmAffinities.forEach((hm, affinityFiles) -> {
                OutputFile hmInput = new OutputFile(groupInput, hm);
                OutputFile hmOutput = addOutput(groupOutput, hm);

                Collection<InputFile> inputs =
                        affinityFiles.stream().map(affinityFile -> addInput(hmInput, affinityFile)).toList();
                bridge.put(hmOutput, inputs);
            });
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            final TreeSet<RegionWithPayload<Set<String>>> regions;
            try (var reader = new BufferedReader(new FileReader(coOccurringRegions))) {
                regions = reader.lines().map(line -> line.split("\t")).map(
                        split -> new RegionWithPayload<Set<String>>(split[0], Integer.parseInt(split[1]),
                                Integer.parseInt(split[2]),
                                new HashSet<>(Arrays.asList(split[3].split("\\|"))))).filter(
                        region -> region.getPayload().size() > 1).collect(Collectors.toCollection(TreeSet::new));
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }

            final Set<String> transcriptionFactors =
                    regions.stream().flatMap(region -> region.getPayload().stream()).collect(Collectors.toSet());

            bridge.forEach((output, inputs) -> add(() -> {
                TreeSet<RegionWithPayload<Map<String, List<Double>>>> regionTfAffinities =
                        inputs.stream().flatMap(input -> {
                            try {
                                List<String> lines = readLines(input);
                                return lines.stream().filter(line -> !line.isBlank()).map(
                                        line -> line.split("\t")).filter(
                                        split -> transcriptionFactors.contains(split[0])).map(split -> {
                                    String region = split[6].substring(1);
                                    double affinity = Double.parseDouble(split[1]);

                                    String[] chromSplit = region.split(":");
                                    String chromosome = chromSplit[0];

                                    String[] posSplit = chromSplit[1].split("-");
                                    int start = Integer.parseInt(posSplit[0]);
                                    int end = Integer.parseInt(posSplit[1]);

                                    RegionWithPayload<Set<String>> searchRegion =
                                            new RegionWithPayload<>(chromosome, start, end, new HashSet<>());

                                    RegionWithPayload<Set<String>> lower = regions.floor(searchRegion);
                                    RegionWithPayload<Set<String>> higher = regions.ceiling(searchRegion);

                                    RegionWithPayload<Set<String>> overlap =
                                            (lower != null && lower.overlaps(searchRegion)) ? lower :
                                                    (higher != null && higher.overlaps(searchRegion)) ? higher : null;

                                    if (overlap == null) {
                                        return null;
                                    }

                                    return Pair.of(overlap, Pair.of(split[0], affinity));
                                });
                            } catch (IOException e) {
                                throw new UncheckedIOException(e);
                            }
                        }).filter(Objects::nonNull).collect(Collectors.groupingBy(Pair::getKey,
                                Collectors.mapping(Pair::getValue, Collectors.toList()))).entrySet().stream().map(
                                entry -> {
                                    RegionWithPayload<Set<String>> region = entry.getKey();
                                    Map<String, List<Double>> tfAffinities = entry.getValue().stream().collect(
                                            Collectors.groupingBy(Pair::getKey,
                                                    Collectors.mapping(Pair::getValue, Collectors.toList())));

                                    return new RegionWithPayload<>(region.getChromosome(), region.getStart(),
                                            region.getEnd(), tfAffinities);
                                }).collect(Collectors.toCollection(TreeSet::new));


                Map<String, BufferedWriter> tfWriters = regionTfAffinities.stream().flatMap(
                        region -> region.getPayload().keySet().stream()).distinct().collect(
                        Collectors.toMap(tf -> tf, tf -> {
                            try {
                                File file = new File(output, tf + ".tsv");
                                makeSureFileExists(file);
                                return new BufferedWriter(new FileWriter(file));
                            } catch (IOException e) {
                                throw new UncheckedIOException(e);
                            }
                        }));

                int i = 0;
                for (RegionWithPayload<Map<String, List<Double>>> region : regionTfAffinities) {
                    Map<String, List<Double>> tfAffinities = region.getPayload();
                    for (String tf : tfAffinities.keySet()) {
                        try {
                            OptionalDouble affinity =
                                    tfAffinities.get(tf).stream().mapToDouble(Double::doubleValue).average();
                            if (affinity.isPresent()) {
                                tfWriters.get(tf).write(i + "\t" + affinity.getAsDouble() + "\n");
                            }
                        } catch (IOException e) {
                            throw new UncheckedIOException(e);
                        }
                    }
                    i++;
                }

                tfWriters.values().forEach(writer -> {
                    try {
                        writer.close();
                    } catch (IOException e) {
                        throw new UncheckedIOException(e);
                    }
                });

                return true;
            }));
        }};
    }
}
