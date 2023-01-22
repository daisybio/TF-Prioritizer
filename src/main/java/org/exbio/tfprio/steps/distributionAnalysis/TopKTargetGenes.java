package org.exbio.tfprio.steps.distributionAnalysis;

import org.apache.commons.lang3.tuple.Pair;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import static org.exbio.pipejar.util.FileManagement.readLines;
import static org.exbio.tfprio.steps.plots.TopKTargetGenes.getTfTopGeneAffinities;
import static org.exbio.tfprio.steps.plots.TopKTargetGenes.writeAffinities;

public class TopKTargetGenes extends ExecutableStep<Configs> {
    public final Map<String, Map<String, OutputFile>> outputFiles = new HashMap<>();
    private final InputFile dcgScoreFile;
    private final Map<InputFile, OutputFile> bridge = new HashMap<>();
    private final RequiredConfig<Integer> topKTargetGenes = new RequiredConfig<>(configs.plots.topKTargetGenes);


    public TopKTargetGenes(Configs configs, OutputFile dcgScoreFile,
                           Map<String, Map<String, OutputFile>> groupHmMeanAffinities) {
        super(configs, false, groupHmMeanAffinities.values().stream().flatMap(sub -> sub.values().stream()).toList(),
                dcgScoreFile);

        this.dcgScoreFile = addInput(dcgScoreFile);

        groupHmMeanAffinities.forEach((group, hmMeanAffinities) -> {
            OutputFile inGroup = new OutputFile(inputDirectory, group);
            OutputFile outGroup = new OutputFile(outputDirectory, group);

            hmMeanAffinities.forEach((hm, affinityFile) -> {
                InputFile inHm = addInput(inGroup, affinityFile);
                OutputFile outHm = addOutput(outGroup, hm);

                outputFiles.computeIfAbsent(group, x -> new HashMap<>()).put(hm, outHm);
                this.bridge.put(inHm, outHm);
            });
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            Set<String> tfNames;
            try {
                tfNames = readLines(dcgScoreFile).stream().skip(1).map(line -> line.split("\t")[0]).collect(
                        Collectors.toSet());
            } catch (IOException e) {
                throw new RuntimeException(e);
            }

            bridge.forEach((inHm, outHm) -> add(() -> {
                Map<String, List<Pair<String, Double>>> tfTopGeneAffinity =
                        getTfTopGeneAffinities(inHm, tfNames, topKTargetGenes.get());

                tfNames.forEach(tf -> {
                    File tfFile = new File(outHm, tf + ".tsv");

                    writeAffinities(tfFile, tfTopGeneAffinity.get(tf));
                });

                return true;
            }));
        }};
    }
}
