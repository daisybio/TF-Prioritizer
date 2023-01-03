package org.exbio.tfprio.steps.distributionAnalysis;

import org.apache.commons.lang3.tuple.Pair;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.pipejar.util.FileFilters.Filters;

import java.io.File;
import java.util.*;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.FileManagement.softLink;

public class ExtractStats extends ExecutableStep {
    public final Map<String, OutputFile> outputFiles = new HashMap<>();

    private final Map<File, OutputFile> bridge = new HashMap<>();

    public ExtractStats(Map<String, OutputFile> plotOutputFiles) {
        super(false, plotOutputFiles.get("stats"));

        InputFile statsDirectory = addInput(plotOutputFiles.get("stats"));
        Arrays.stream(Objects.requireNonNull(statsDirectory.listFiles(Filters.fileFilter))).map(
                file -> Pair.of(file, addOutput(file.getName()))).forEach(pair -> {
            bridge.put(pair.getLeft(), pair.getRight());
            outputFiles.put(pair.getLeft().getName(), pair.getRight());
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((input, output) -> add(() -> {
                softLink(output, input);

                return true;
            }));
        }};
    }
}
