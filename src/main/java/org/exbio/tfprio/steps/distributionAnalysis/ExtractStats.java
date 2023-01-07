package org.exbio.tfprio.steps.distributionAnalysis;

import org.apache.commons.lang3.tuple.Pair;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;

import java.io.File;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.FileManagement.softLink;

public class ExtractStats extends ExecutableStep {
    public final Map<String, OutputFile> outputFiles = new HashMap<>();

    private final Map<File, OutputFile> bridge = new HashMap<>();

    public ExtractStats(Map<String, Pair<OutputFile, OutputFile>> plotOutputFiles) {
        super(false, plotOutputFiles.values().stream().map(Pair::getRight).toList());

        plotOutputFiles.forEach((hm, output) -> {
            OutputFile outputStats = output.getRight();
            OutputFile outputHm = addOutput(hm);
            InputFile inputHm = addInput(outputStats);
            outputFiles.put(hm, outputHm);
            bridge.put(inputHm, outputHm);
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
