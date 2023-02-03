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
import java.util.stream.IntStream;

import static org.exbio.pipejar.util.FileManagement.makeSureFileExists;

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
            final List<RegionWithPayload<Set<String>>> regions;
            try (var reader = new BufferedReader(new FileReader(coOccurringRegions))) {
                regions = reader.lines().map(line -> line.split("\t")).map(
                        split -> new RegionWithPayload<Set<String>>(split[0], Integer.parseInt(split[1]),
                                Integer.parseInt(split[2]),
                                new HashSet<>(Arrays.asList(split[3].split("\\|"))))).filter(
                        region -> region.getPayload().size() > 1).sorted().toList();
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }

            final Set<String> transcriptionFactors =
                    regions.stream().flatMap(region -> region.getPayload().stream()).collect(Collectors.toSet());

            bridge.forEach((output, inputs) -> add(() -> {
                logger.trace("Loading affinities");
                Map<String, TreeSet<RegionWithPayload<Double>>> tfAffinities = inputs.stream().flatMap(input -> {
                    try (var reader = new BufferedReader(new FileReader(input))) {
                        return reader.lines().filter(line -> !line.isBlank()).map(line -> line.split("\t")).filter(
                                split -> transcriptionFactors.contains(split[0])).map(split -> {
                            logger.trace("Parsing line: " + Arrays.toString(split));
                            String region = split[6].substring(1);
                            double affinity = Double.parseDouble(split[1]);

                            String[] chromSplit = region.split(":");
                            String chromosome = chromSplit[0];

                            String[] posSplit = chromSplit[1].split("-");
                            int start = Integer.parseInt(posSplit[0]);
                            int end = Integer.parseInt(posSplit[1]);

                            logger.trace("Parsed line: " + chromosome + ":" + start + "-" + end + " " + affinity);

                            return Pair.of(split[0], new RegionWithPayload<>(chromosome, start, end, affinity));
                        });
                    } catch (IOException e) {
                        logger.warn("Failed to read file: " + input);
                        throw new UncheckedIOException(e);
                    }
                }).collect(Collectors.groupingBy(Pair::getKey,
                        Collectors.mapping(Pair::getValue, Collectors.toCollection(TreeSet::new))));

                logger.trace("Opening writers");
                Map<String, BufferedWriter> tfWriters =
                        tfAffinities.keySet().stream().collect(Collectors.toMap(tf -> tf, tf -> {
                            try {
                                File file = new File(output, tf + ".tsv");
                                makeSureFileExists(file);
                                return new BufferedWriter(new FileWriter(file));
                            } catch (IOException e) {
                                logger.warn("Failed to create file: " + output + "/" + tf + ".tsv");
                                throw new UncheckedIOException(e);
                            }
                        }));

                logger.trace("Writing");
                IntStream.range(0, regions.size()).forEach(i -> {
                    RegionWithPayload<Set<String>> region = regions.get(i);

                    RegionWithPayload<Double> searchRegion =
                            new RegionWithPayload<>(region.getChromosome(), region.getStart(), region.getEnd(), null);

                    Set<String> regionTfs = region.getPayload();
                    regionTfs.forEach(tf -> {
                        TreeSet<RegionWithPayload<Double>> affinities = tfAffinities.get(tf);

                        double averageAffinity = affinities.subSet(searchRegion, searchRegion).stream().mapToDouble(
                                RegionWithPayload::getPayload).average().orElse(0.0);

                        try {
                            tfWriters.get(tf).write(i + "\t" + averageAffinity + "\n");
                        } catch (IOException e) {
                            logger.warn("Failed to write to file: " + output + "/" + tf + ".tsv");
                            throw new UncheckedIOException(e);
                        }
                    });
                });

                logger.trace("Closing writers");
                tfWriters.values().forEach(writer -> {
                    try {
                        writer.close();
                    } catch (IOException e) {
                        logger.warn("Failed to close file");
                        throw new UncheckedIOException(e);
                    }
                });

                return true;
            }));
        }};
    }
}
