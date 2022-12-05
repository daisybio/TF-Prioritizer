package org.exbio.tfprio.steps.tGene;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import static org.exbio.pipejar.util.FileManagement.deleteFileStructure;
import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class TGene extends ExecutableStep {
    public final Map<String, Map<String, Collection<OutputFile>>> outputFiles = new HashMap<>();
    private final Map<InputFile, OutputFile> bridge = new HashMap<>();
    private final InputFile gtfFile;
    private final RequiredConfig<Boolean> noClosestLocus = new RequiredConfig<>(Configs.tGene.noClosestLocus);
    private final RequiredConfig<Boolean> noClosestTss = new RequiredConfig<>(Configs.tGene.noClosestTss);
    private final RequiredConfig<Integer> maxLinkDistance = new RequiredConfig<>(Configs.tGene.maxLinkDistance);
    private final RequiredConfig<Double> pValue = new RequiredConfig<>(Configs.tGene.pValue);

    private final RequiredConfig<File> executable = new RequiredConfig<>(Configs.tGene.executable);

    public TGene(Map<String, Map<String, Collection<OutputFile>>> peakFiles, OutputFile gtfFile) {
        super(false, peakFiles.values().stream().flatMap(
                stringCollectionMap -> stringCollectionMap.values().stream()).flatMap(Collection::stream).collect(
                Collectors.toSet()), gtfFile);

        this.gtfFile = addInput(gtfFile);

        peakFiles.forEach((group, hmMap) -> {
            outputFiles.put(group, new HashMap<>());
            OutputFile d_groupOut = addOutput(group);
            hmMap.forEach((hm, sampleFiles) -> {
                outputFiles.get(group).put(hm, new HashSet<>());
                OutputFile d_hmOut = addOutput(d_groupOut, hm);

                sampleFiles.forEach(sampleFile -> {
                    InputFile inputFile = addInput(sampleFile);
                    final OutputFile sampleDirectory;
                    if (sampleFiles.size() > 1) {
                        sampleDirectory = new OutputFile(d_hmOut,
                                sampleFile.getName().substring(0, sampleFile.getName().lastIndexOf('.')));
                    } else {
                        sampleDirectory = d_hmOut;
                    }
                    OutputFile links = addOutput(sampleDirectory, sampleDirectory.getName() + ".tsv");
                    outputFiles.get(group).get(hm).add(links);
                    bridge.put(inputFile, links);
                });
            });
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((inputFile, outputFile) -> add(() -> {
                String command = executable.get().getAbsolutePath() + " " + inputFile.getAbsolutePath() + " " +
                        gtfFile.getAbsolutePath() + " -oc " + outputFile.getParentFile().getAbsolutePath();

                if (noClosestLocus.get()) {
                    command += " --no-closest-locus";
                }
                if (noClosestTss.get()) {
                    command += " --no-closest-tss";
                }

                command += " --max-link-distances " + maxLinkDistance.get();
                command += " --max-pvalue " + pValue.get();

                executeAndWait(command, true);

                File originalOutput = new File(outputFile.getParentFile(), "links.tsv");

                try (BufferedReader reader = new BufferedReader(new java.io.FileReader(originalOutput));
                     BufferedWriter writer = new BufferedWriter(new java.io.FileWriter(outputFile))) {
                    String line;
                    while ((line = reader.readLine()) != null) {
                        writer.write(line);
                        writer.newLine();
                    }
                }

                deleteFileStructure(originalOutput);

                return true;
            }));
        }};
    }
}
