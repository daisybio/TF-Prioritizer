package org.exbio.tfprio.steps.peakFiles;

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
import java.util.stream.Stream;

import static org.exbio.pipejar.util.FileManagement.readLines;

/**
 *
 */
public class MixMutuallyExclusive extends ExecutableStep<Configs> {
    public final Map<String, Map<String, Collection<OutputFile>>> outputFiles = new HashMap<>();
    private final Map<String, Map<String, Collection<InputFile>>> inputFiles = new HashMap<>();
    private final RequiredConfig<Boolean> differentialPeakSignals =
            new RequiredConfig<>(configs.mixOptions.differentialPeakSignals);

    public MixMutuallyExclusive(Configs configs, Map<String, Map<String, Collection<OutputFile>>> dependencies) {
        super(configs, false, dependencies.values().stream().flatMap(
                stringCollectionMap -> stringCollectionMap.values().stream()).flatMap(Collection::stream).collect(
                Collectors.toSet()));
        dependencies.forEach((group, hmMap) -> {
            outputFiles.put(group, new HashMap<>());
            inputFiles.put(group, new HashMap<>());
            OutputFile d_groupOut = addOutput(group);
            hmMap.forEach((hm, sampleFiles) -> {
                if (sampleFiles.size() > 1) {
                    throw new IllegalArgumentException("MixMutuallyExclusive can only handle one sample per hm");
                }

                outputFiles.get(group).put(hm, new HashSet<>());
                inputFiles.get(group).put(hm, new HashSet<>());
                OutputFile d_hmOut = addOutput(d_groupOut, hm);

                sampleFiles.forEach(sampleFile -> {
                    InputFile inputFile = addInput(sampleFile);
                    OutputFile outputFile = addOutput(d_hmOut, inputFile.getName());
                    outputFiles.get(group).get(hm).add(outputFile);
                    inputFiles.get(group).get(hm).add(inputFile);
                });
            });
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            Map<String, Map<String, List<BroadPeakRegion>>> regions = new HashMap<>();

            inputFiles.forEach((group, hmMap) -> {
                regions.put(group, new HashMap<>());
                hmMap.forEach((hm, sampleFiles) -> {
                    regions.get(group).put(hm, new ArrayList<>());
                    sampleFiles.forEach(inputFile -> {
                        try {
                            regions.get(group).get(hm).addAll(
                                    readLines(inputFile).stream().map(BroadPeakRegion::new).toList());
                        } catch (IOException e) {
                            throw new RuntimeException(e);
                        }
                    });
                });
            });

            outputFiles.forEach((group, hmMap) -> hmMap.forEach((hm, sampleFiles) -> sampleFiles.forEach(outputFile ->
                    // Here the callables are added
                    add(() -> {
                        TreeSet<BroadPeakRegion> others = regions.entrySet().stream()
                                                                 // Investigate all groups
                                                                 .flatMap(groupEntry ->
                                                                         // Investigate all hms
                                                                         groupEntry.getValue().entrySet().stream()
                                                                                   // Keep only the entries, which are not fully equal to the active entry
                                                                                   .filter(hmEntry -> !(
                                                                                           groupEntry.getKey().equals(
                                                                                                   group) &&
                                                                                                   hmEntry.getKey().equals(
                                                                                                           hm))))
                                                                 // Map the remaining entries to their regions
                                                                 .map(Map.Entry::getValue)
                                                                 // Flatten the stream of lists to a stream of regions
                                                                 .flatMap(Collection::stream).collect(
                                        // Collect to TreeSet
                                        Collectors.toCollection(TreeSet::new));
                        double averageScore =
                                others.stream().mapToDouble(BroadPeakRegion::getScore).average().orElse(0);

                        logger.debug("Start filtering");
                        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
                            regions.get(group).get(hm).stream()
                                   // Filter out all regions, which are not mutually exclusive
                                   .filter(region -> {
                                       // Check if any other region overlaps with the current region
                                       Optional<BroadPeakRegion> match =
                                               // Get elements next to the current region
                                               Stream.of(others.floor(region), others.ceiling(region))
                                                     // Check if they exist and overlap
                                                     .filter(other -> other != null && other.overlaps(region))
                                                     // Overlapping means that this is our match
                                                     .findFirst();

                                       if (match.isEmpty()) {
                                           // If there is no match, the region is mutually exclusive
                                           return true;
                                       } else if (differentialPeakSignals.get()) {
                                           // If this option is enabled, scores are considered
                                           double scoreDifference =
                                                   Math.abs(region.getScore() - match.get().getScore());

                                           return scoreDifference > averageScore;
                                       } else {
                                           // If none of the other conditions are met, the region is not mutually exclusive
                                           return false;
                                       }
                                   })
                                   // Sort and save the regions
                                   .sorted().forEachOrdered(region -> {
                                       try {
                                           writer.write(region.toString());
                                           writer.newLine();
                                       } catch (IOException e) {
                                           throw new RuntimeException(e);
                                       }
                                   });
                        }
                        return true;
                    }))));
        }};
    }
}
