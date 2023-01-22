package org.exbio.tfprio.steps.peakFiles;

import org.apache.commons.lang3.tuple.ImmutablePair;
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

/**
 * This step takes the output of the {@link CheckChromosomes} step and
 * combines the files into one file per group and histone mark. Depending on the configuration,
 * either all overlapping peaks are merged or only the ones that occur more than
 * {@link org.exbio.tfprio.configs.Configs#mixOptions#minOccurrence} times.
 */
public class MixSamples extends ExecutableStep<Configs> {
    public final Map<String, Map<String, Collection<OutputFile>>> outputFiles = new HashMap<>();
    private final RequiredConfig<Integer> minOccurrence = new RequiredConfig<>(configs.mixOptions.minOccurrence);
    private final Map<OutputFile, Collection<InputFile>> bridge = new HashMap<>();

    public MixSamples(Configs configs, Map<String, Map<String, Collection<OutputFile>>> dependencies) {
        super(configs, false, dependencies.values().stream().flatMap(
                stringCollectionMap -> stringCollectionMap.values().stream()).flatMap(Collection::stream).collect(
                Collectors.toSet()));
        dependencies.forEach((group, hmMap) -> {
            outputFiles.put(group, new HashMap<>());
            OutputFile d_groupOut = addOutput(group);
            hmMap.forEach((hm, sampleFiles) -> {
                OutputFile f_output = addOutput(d_groupOut, group + "_" + hm + ".broadPeak");
                outputFiles.get(group).put(hm, new HashSet<>() {{
                    add(f_output);
                }});
                bridge.put(f_output, new HashSet<>());
                sampleFiles.forEach(sampleFile -> {
                    InputFile inputFile = addInput(sampleFile);
                    bridge.get(f_output).add(inputFile);
                });
            });
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((outputFile, inputFiles) -> add(() -> {
                // Collect all regions from all input files
                Map<InputFile, List<Region>> sortedRegions = new HashMap<>();
                inputFiles.forEach(inputFile -> {
                    try {
                        sortedRegions.put(inputFile, readLines(inputFile).stream().map(Region::new).sorted().toList());
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                });

                Map<InputFile, Integer> currentIndices = new HashMap<>() {{
                    inputFiles.forEach(inputFile -> put(inputFile, 0));
                }};

                SortedMap<Region, Collection<Region>> overlappingRegions = new TreeMap<>();

                Optional<Region> optionalRegion;

                while (
                    // Find next region to process
                        (optionalRegion = sortedRegions.entrySet().stream()
                                                       // Only investigate map entries with remaining regions
                                                       .filter(entry -> entry.getValue().size() >
                                                               currentIndices.get(entry.getKey()))
                                                       // Map to InputFile-Region pairs for current index
                                                       .map(entry -> new ImmutablePair<>(entry.getKey(),
                                                               entry.getValue().get(
                                                                       currentIndices.get(entry.getKey()))))
                                                       // Find minimum over all files
                                                       .min(Comparator.comparing(pair -> pair.right))
                                                       // Increase index for file with minimum and map to region
                                                       .map(pair -> {
                                                           currentIndices.put(pair.left,
                                                                   currentIndices.get(pair.left) + 1);
                                                           return pair.right;
                                                       })).isPresent()) {

                    Region region = optionalRegion.get();

                    // Only entered in first iteration
                    if (overlappingRegions.isEmpty()) {
                        overlappingRegions.put(new Region(region.getChromosome(), region.getStart(), region.getEnd()),
                                new HashSet<>());
                    }

                    // Check out if this region overlaps with last investigated region
                    Region lastRegion = overlappingRegions.lastKey();
                    if (lastRegion.overlaps(region)) {
                        // Overlap: Add region to overlapping group
                        overlappingRegions.get(lastRegion).add(region);
                        lastRegion.merge(region);
                    } else {
                        // No overlap: Create new group from this region
                        overlappingRegions.put(new Region(region.getChromosome(), region.getStart(), region.getEnd()),
                                new HashSet<>() {{
                                    add(region);
                                }});
                    }
                }

                final int minCount =
                        minOccurrence.isSet() ? Math.min(minOccurrence.get(), inputFiles.size()) : inputFiles.size();

                List<Region> mergedRegions = overlappingRegions.values().stream()
                                                               // Keep only groups with at least minCount regions
                                                               .filter(group -> group.size() >= minCount)
                                                               // Merge remaining into one region per overlapping group
                                                               .map(Region::mergeMulti)
                                                               // Get sorted list
                                                               .sorted().toList();

                // Write the output file
                try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
                    for (Region region : mergedRegions) {
                        writer.write(region.toString());
                        writer.newLine();
                    }
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }

                return true;
            }));
        }};
    }
}
