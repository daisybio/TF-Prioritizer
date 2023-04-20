package org.exbio.tfprio.steps.ehmm;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.File;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import static org.exbio.pipejar.util.FileManagement.softLink;

public class InitBamFiles extends ExecutableStep<Configs> {
    public final Map<OutputFile, Map<OutputFile, Collection<OutputFile>>> outputFiles = new HashMap<>();
    private final RequiredConfig<File> inputBamDir = new RequiredConfig<>(configs.ehmm.bamDirectory);
    private final Map<File, OutputFile> bridge = new HashMap<>();
    
    public InitBamFiles(Configs configs) {
        super(configs);
        InputFile inputBamDirectory = addInput(inputBamDir);
        Arrays.stream(Objects.requireNonNull(inputBamDirectory.listFiles())).filter(File::isDirectory).forEach(d_group -> {
            OutputFile d_groupOut = new OutputFile(outputDirectory, d_group.getName());
            Arrays.stream(Objects.requireNonNull(d_group.listFiles())).filter(File::isDirectory).forEach(d_hm -> {
                OutputFile d_hmOut = new OutputFile(d_groupOut, d_hm.getName());
                List<File> bamFiles = Arrays.stream(Objects.requireNonNull(d_hm.listFiles()))
                        .filter(File::exists).filter(f -> f.getAbsolutePath().endsWith(".bam"))
                        .collect(Collectors.toList());
                bamFiles.forEach(f_sample -> {
                    OutputFile f_sampleOut = addOutput(d_hmOut, f_sample.getName());
                    bridge.put(f_sample, f_sampleOut);
                    outputFiles.computeIfAbsent(d_groupOut, k -> new HashMap<>()).computeIfAbsent(d_hmOut,
                            k -> new HashSet<>()).add(f_sampleOut);
                });
            });
        });
    }

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
