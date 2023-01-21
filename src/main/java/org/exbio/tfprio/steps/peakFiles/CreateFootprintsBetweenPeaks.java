package org.exbio.tfprio.steps.peakFiles;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;
import org.exbio.tfprio.lib.Region;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import static org.exbio.pipejar.util.FileManagement.readLines;

public class CreateFootprintsBetweenPeaks extends ExecutableStep<Configs> {
    private final Map<InputFile, OutputFile> bridge = new HashMap<>();
    private final RequiredConfig<String> tfBindingSiteSearch =
            new RequiredConfig<>(configs.mixOptions.tfBindingSiteSearch);
    private final RequiredConfig<Integer> maxDistance = new RequiredConfig<>(configs.mixOptions.maxSpaceBetweenPeaks);
    public Map<String, Map<String, Collection<OutputFile>>> outputFiles = new HashMap<>();

    public CreateFootprintsBetweenPeaks(Configs configs,
                                        Map<String, Map<String, Collection<OutputFile>>> dependencies) {
        super(configs, false, dependencies.values().stream().flatMap(
                stringCollectionMap -> stringCollectionMap.values().stream()).flatMap(Collection::stream).collect(
                Collectors.toSet()));

        dependencies.forEach((group, hmMap) -> {
            outputFiles.put(group, new HashMap<>());
            OutputFile d_groupOut = addOutput(group);
            hmMap.forEach((hm, sampleFiles) -> {
                outputFiles.get(group).put(hm, new HashSet<>());
                OutputFile d_hmOut = addOutput(d_groupOut, hm);

                sampleFiles.forEach(sampleFile -> {
                    InputFile inputFile = addInput(sampleFile);
                    OutputFile outputFile = addOutput(d_hmOut, inputFile.getName());
                    bridge.put(inputFile, outputFile);
                    outputFiles.get(group).get(hm).add(outputFile);
                });
            });
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((inputFile, outputFile) -> add(() -> {
                final List<Region> regions = readLines(inputFile).stream().map(Region::new).toList();

                if (regions.isEmpty()) {
                    logger.warn("No regions found in file " + inputFile);
                    return true;
                }

                final List<Region> processedRegions;

                if (tfBindingSiteSearch.get().equals("BETWEEN")) {
                    List<Set<Region>> groupedRegions = new ArrayList<>();

                    groupedRegions.add(new HashSet<>() {{
                        add(regions.get(0));
                    }});

                    Region previousRegion = regions.get(0);

                    for (Region region : regions.subList(1, regions.size())) {
                        if (region.getStart() - previousRegion.getEnd() <= maxDistance.get()) {
                            groupedRegions.get(groupedRegions.size() - 1).add(region);
                        } else {
                            groupedRegions.add(new HashSet<>() {{
                                add(region);
                            }});
                        }
                        previousRegion = region;
                    }
                    processedRegions = groupedRegions.stream().map(Region::mergeMulti).toList();
                } else {
                    // EXCL_BETWEEN
                    Region previousRegion = regions.get(0);

                    processedRegions = new ArrayList<>();

                    for (Region region : regions.subList(1, regions.size())) {
                        if ((region.getChromosome().equals(previousRegion.getChromosome())) &&
                                (region.getStart() - previousRegion.getEnd() <= maxDistance.get())) {
                            processedRegions.add(Region.between(previousRegion, region));
                        }
                        previousRegion = region;
                    }
                }


                try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
                    processedRegions.forEach(region -> {
                        try {
                            writer.write(region.toString());
                            writer.newLine();
                        } catch (IOException e) {
                            throw new RuntimeException(e);
                        }
                    });
                }

                return true;
            }));
        }};
    }
}
