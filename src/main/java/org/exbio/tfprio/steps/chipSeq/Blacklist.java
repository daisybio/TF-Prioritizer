package org.exbio.tfprio.steps.chipSeq;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;
import org.exbio.tfprio.lib.BedRegion;
import org.exbio.tfprio.lib.BroadPeakRegion;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import static org.exbio.pipejar.util.FileManagement.readLines;

public class Blacklist extends ExecutableStep {
    public final Map<String, Map<String, Collection<OutputFile>>> outputFiles = new HashMap<>();
    private final Map<InputFile, OutputFile> bridge = new HashMap<>();
    private final RequiredConfig<File> blacklist = new RequiredConfig<>(Configs.mixOptions.blackListPath);

    public Blacklist(Map<String, Map<String, Collection<OutputFile>>> dependencies) {
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
            try {
                // Load all blacklisted regions
                Collection<BedRegion> blacklistRegions =
                        readLines(blacklist.get()).stream().map(BedRegion::new).toList();

                bridge.forEach((inputFile, outputFile) -> add(() -> {
                    try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
                        // Read all regions from the input file
                        readLines(inputFile).stream().map(BroadPeakRegion::new)
                                            // Filter out all blacklisted regions
                                            .filter(broadPeakRegion -> blacklistRegions.stream().noneMatch(
                                                    broadPeakRegion::overlaps))
                                            // Get string representations of remaining regions
                                            .map(BroadPeakRegion::toString)
                                            // Write them to the output file
                                            .forEach(line -> {
                                                try {
                                                    writer.write(line);
                                                    writer.newLine();
                                                } catch (IOException e) {
                                                    throw new RuntimeException(e);
                                                }
                                            });
                    }
                    return true;
                }));
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }};
    }
}
