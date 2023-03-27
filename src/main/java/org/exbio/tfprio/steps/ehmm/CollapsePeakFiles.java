package org.exbio.tfprio.steps.ehmm;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;


public class CollapsePeakFiles extends ExecutableStep<Configs> {
    public final Map<String, OutputFile> outputFiles = new HashMap<>();
    private final Map<OutputFile, Collection<InputFile>> bridge = new HashMap<>();

    public CollapsePeakFiles(Configs configs, Map<String, Map<String, Collection<OutputFile>>> peakFiles) {
        super(configs, false, peakFiles.values().stream().flatMap(
                stringCollectionMap -> stringCollectionMap.values().stream()).flatMap(Collection::stream).collect(
                Collectors.toSet()));
        peakFiles.forEach((group, hmMap) -> {
            OutputFile dir = addOutput(group);
            OutputFile out = addOutput(dir, "collapsedPeaks.bed");
            outputFiles.put(group, out);
            // ignore hm groups and collapse all bed files to a list
            Collection<InputFile> files = hmMap.values().stream().
                    flatMap(Collection::stream)
                    .map(this::addInput)
                    .collect(Collectors.toSet());
            bridge.put(out, files);
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                bridge.forEach((out, peakFiles) -> peakFiles
                        .forEach(file -> {
                            try {
                                Files.write(out.toPath(),
                                    (Iterable<String>) Files.lines(file.toPath())::iterator,
                                    StandardOpenOption.APPEND);
                            } catch (IOException e) {
                                throw new RuntimeException("Failed to collapse peaks when adding: "
                                        + file.getName() + " in group: " + out.getName());
                            }
                        }));
                return true;
            });
        }};
    }
}
