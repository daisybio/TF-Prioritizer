package org.exbio.tfprio.steps.TEPIC;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;
import org.exbio.tfprio.lib.Region;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

public class TEPICPostprocessing extends ExecutableStep {
    public final Map<String, Map<String, Collection<OutputFile>>> outputFiles = new HashMap<>();
    private final Map<InputFile, OutputFile> bridge = new HashMap<>();
    private final RequiredConfig<Double> affinityCutoff = new RequiredConfig<>(Configs.tepic.affinityCutoff);
    private final RequiredConfig<String> sequenceFileName = new RequiredConfig<>(Configs.tepic.sequenceFileName);

    public TEPICPostprocessing(Map<String, Map<String, Collection<OutputFile>>> tepicResults) {
        super(false,
                tepicResults.values().stream().flatMap(x -> x.values().stream()).flatMap(Collection::stream).collect(
                        Collectors.toList()));
        tepicResults.forEach((group, hmMap) -> {
            outputFiles.put(group, new HashMap<>());
            OutputFile outputGroup = addOutput(group);
            hmMap.forEach((hm, sampleDirectories) -> {
                outputFiles.get(group).put(hm, new HashSet<>(sampleDirectories.size()));
                OutputFile outputHm = addOutput(outputGroup, hm);
                sampleDirectories.forEach(sampleDirectory -> {
                    InputFile inputDirectory = addInput(sampleDirectory);
                    OutputFile outputDirectory = addOutput(outputHm, sampleDirectory.getName());
                    bridge.put(inputDirectory, outputDirectory);
                    outputFiles.get(group).get(hm).add(outputDirectory);
                });
            });
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((inputDirectory, outputDirectory) -> {
                File sequenceFile = new File(inputDirectory, sequenceFileName.get());

                Map<String, Collection<Region>> tfRegions = new HashMap<>();

                try (BufferedReader reader = new BufferedReader(new FileReader(sequenceFile))) {
                    for (String line = reader.readLine(); line != null; line = reader.readLine()) {
                        if (line.startsWith("#")) {
                            continue;
                        }

                        String[] split = line.split("\t");

                        String tf = split[0];
                        double affinity = Double.parseDouble(split[1]);
                        int start = Integer.parseInt(split[3]);
                        int length = Integer.parseInt(split[4]);
                        String sequence = split[5];

                        if (affinity < affinityCutoff.get()) {
                            continue;
                        }
                    }
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            });
        }};
    }
}
