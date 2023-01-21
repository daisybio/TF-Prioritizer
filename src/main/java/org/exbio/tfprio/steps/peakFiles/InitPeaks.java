package org.exbio.tfprio.steps.peakFiles;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.File;
import java.util.*;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.FileManagement.softLink;

public class InitPeaks extends ExecutableStep<Configs> {
    public final Map<String, Map<String, Collection<OutputFile>>> outputFiles = new HashMap<>();
    private final RequiredConfig<File> inputDirectory = new RequiredConfig<>(configs.inputConfigs.peaks);
    private final Map<File, OutputFile> bridge = new HashMap<>();

    public InitPeaks(Configs configs) {
        super(configs);

        // Register all input files and their corresponding output files
        InputFile input = addInput(inputDirectory);

        Arrays.stream(Objects.requireNonNull(input.listFiles())).forEach(d_group -> {
            String group = d_group.getName();
            OutputFile d_groupOut = new OutputFile(outputDirectory, d_group.getName());

            Arrays.stream(Objects.requireNonNull(d_group.listFiles())).forEach(d_hm -> {
                String hm = d_hm.getName();
                OutputFile d_hmOut = new OutputFile(d_groupOut, d_hm.getName());

                Arrays.stream(Objects.requireNonNull(d_hm.listFiles())).forEach(f_sample -> {
                    OutputFile f_sampleOut = addOutput(d_hmOut, f_sample.getName());
                    bridge.put(f_sample, f_sampleOut);
                    outputFiles.computeIfAbsent(group, k -> new HashMap<>()).computeIfAbsent(hm,
                            k -> new HashSet<>()).add(f_sampleOut);
                });
            });
        });
    }

    //TODO: check if peak files have correct bed format
    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((inputFile, outputFile) -> add(() -> {
                softLink(outputFile, inputFile);
                return true;
            }));
        }};
    }
}