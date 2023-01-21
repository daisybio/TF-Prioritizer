package org.exbio.tfprio.steps.distributionAnalysis;

import org.apache.commons.lang3.tuple.Pair;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.File;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class CreateHeatmaps extends ExecutableStep<Configs> {
    public final Map<String, Map<String, OutputFile>> outputFiles = new HashMap<>();
    private final InputFile script;
    private final Map<OutputFile, Pair<InputFile, Pair<InputFile, InputFile>>> bridge = new HashMap<>();
    private final InputFile batchFile;
    private final InputFile ensgSymbol;

    public CreateHeatmaps(Configs configs, Map<String, OutputFile> normalizedCounts,
                          Map<String, Map<String, OutputFile>> targetGenes, OutputFile batchFile,
                          OutputFile ensgSymbol) {
        super(configs, false, normalizedCounts.values(),
                targetGenes.values().stream().flatMap(map -> map.values().stream()).toList(), List.of(batchFile));

        script = addInput(getClass().getResourceAsStream("heatmaps.R"), "heatmaps.R");
        this.batchFile = addInput(batchFile);
        this.ensgSymbol = addInput(ensgSymbol);

        Collection<String> hms =
                targetGenes.values().stream().map(Map::keySet).flatMap(Collection::stream).collect(Collectors.toSet());

        Map<String, Map<String, InputFile>> inputTargetGenes =
                targetGenes.entrySet().stream().collect(Collectors.toMap(Map.Entry::getKey, entry -> {
                    OutputFile outputGroup = new OutputFile(inputDirectory, entry.getKey());
                    return entry.getValue().entrySet().stream().collect(Collectors.toMap(Map.Entry::getKey,
                            e -> addInput(new OutputFile(outputGroup, e.getKey()), e.getValue())));
                }));

        normalizedCounts.forEach((pairing, normalizedCount) -> {
            InputFile inputCounts = addInput(normalizedCount);

            String[] groups = pairing.split("_");
            String group1 = groups[0];
            String group2 = groups[1];

            OutputFile outPairing = new OutputFile(outputDirectory, pairing);

            hms.forEach(hm -> {
                InputFile targetGenes1 = inputTargetGenes.get(group1).get(hm);
                InputFile targetGenes2 = inputTargetGenes.get(group2).get(hm);

                if (targetGenes1 == null || targetGenes2 == null) {
                    return;
                }

                OutputFile output = addOutput(outPairing, hm);

                bridge.put(output, Pair.of(inputCounts, Pair.of(targetGenes1, targetGenes2)));
                outputFiles.computeIfAbsent(pairing, k -> new HashMap<>()).put(hm, output);
            });
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((outputDirectory, input) -> {
                InputFile counts = input.getLeft();
                InputFile targetGenesDirectory1 = input.getRight().getLeft();
                InputFile targetGenesDirectory2 = input.getRight().getRight();

                Map<String, Pair<File, File>> tfTargetGenes =
                        Arrays.stream(Objects.requireNonNull(targetGenesDirectory1.listFiles())).collect(
                                Collectors.toMap(f -> f.getName().substring(0, f.getName().lastIndexOf(".")),
                                        file -> Pair.of(file, new File(targetGenesDirectory2, file.getName()))));

                tfTargetGenes.forEach((tf, targetGenesPair) -> add(() -> {
                    File targetGenes1 = targetGenesPair.getLeft();
                    File targetGenes2 = targetGenesPair.getRight();

                    if (!targetGenes1.exists() || !targetGenes2.exists()) {
                        logger.warn("Target genes for " + tf + " not found in " + targetGenesDirectory1 + " or " +
                                targetGenesDirectory2);
                        return false;
                    }

                    File output = new File(outputDirectory, tf + ".png");
                    File outCounts = new File(outputDirectory, tf + ".tsv");

                    String command =
                            "Rscript " + script + " --counts " + counts + " --tg1 " + targetGenes1 + " --tg2 " +
                                    targetGenes2 + " --heatmap " + output + " --outCounts " + outCounts + " --groups " +
                                    batchFile + " --ensgSymbol " + ensgSymbol;

                    executeAndWait(command, false);

                    return true;
                }));
            });
        }};
    }
}
