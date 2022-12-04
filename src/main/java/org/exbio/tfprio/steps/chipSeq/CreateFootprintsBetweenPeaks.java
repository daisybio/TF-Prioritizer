package org.exbio.tfprio.steps.chipSeq;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;
import org.exbio.tfprio.lib.BroadPeakRegion;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import static org.exbio.pipejar.util.FileManagement.readLines;

public class CreateFootprintsBetweenPeaks extends ExecutableStep {
    private final Map<InputFile, OutputFile> bridge = new HashMap<>();
    private final RequiredConfig<String> tfBindingSiteSearch =
            new RequiredConfig<>(Configs.mixOptions.tfBindingSiteSearch);
    private final RequiredConfig<Integer> maxDistance = new RequiredConfig<>(Configs.mixOptions.maxSpaceBetweenPeaks);
    public Map<String, Map<String, Collection<OutputFile>>> outputFiles = new HashMap<>();

    public CreateFootprintsBetweenPeaks(Map<String, Map<String, Collection<OutputFile>>> dependencies) {
        super(false, dependencies.values().stream().flatMap(
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
                final List<BroadPeakRegion> regions = readLines(inputFile).stream().map(BroadPeakRegion::new).toList();
                final List<BroadPeakRegion> processedRegions;

                if (tfBindingSiteSearch.get().equals("BETWEEN")) {
                    List<Set<BroadPeakRegion>> groupedRegions = new ArrayList<>();

                    groupedRegions.add(new HashSet<>() {{
                        add(regions.get(0));
                    }});

                    BroadPeakRegion previousRegion = regions.get(0);

                    for (BroadPeakRegion region : regions.subList(1, regions.size())) {
                        if (region.getStart() - previousRegion.getEnd() <= maxDistance.get()) {
                            groupedRegions.get(groupedRegions.size() - 1).add(region);
                        } else {
                            groupedRegions.add(new HashSet<>() {{
                                add(region);
                            }});
                        }
                        previousRegion = region;
                    }
                    processedRegions = groupedRegions.stream().map(BroadPeakRegion::merge).toList();
                } else {
                    // EXCL_BETWEEN
                    BroadPeakRegion previousRegion = regions.get(0);

                    processedRegions = new ArrayList<>();

                    for (BroadPeakRegion region : regions.subList(1, regions.size())) {
                        if (region.getStart() - previousRegion.getEnd() <= maxDistance.get()) {
                            processedRegions.add(BroadPeakRegion.between(previousRegion, region));
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
                            e.printStackTrace();
                        }
                    });
                }

                return true;
            }));
        }};
    }
}
