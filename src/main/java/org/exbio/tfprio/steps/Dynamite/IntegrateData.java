package org.exbio.tfprio.steps.Dynamite;

import org.apache.commons.lang3.tuple.Pair;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;

import java.util.*;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class IntegrateData extends ExecutableStep {
    public final Map<String, Map<String, OutputFile>> outputFiles = new HashMap<>();
    private final InputFile script;
    private final Map<OutputFile, Pair<InputFile, InputFile>> bridge = new HashMap<>();

    public IntegrateData(Map<String, Map<String, OutputFile>> affinityRatios, Collection<OutputFile> diffExpression) {
        super(false, affinityRatios.values().stream().flatMap(
                stringCollectionMap -> stringCollectionMap.values().stream()).toList(), diffExpression);
        script = addInput(
                getClass().getResourceAsStream("/org/exbio/tfprio/steps/TEPIC/DYNAMITE/Scripts/integrateData.py"),
                "integrate_data.py");

        OutputFile inputDiffExpression = new OutputFile(inputDirectory, "diffExpression");
        OutputFile inputAffinityRatios = new OutputFile(inputDirectory, "affinityRatios");

        Map<String, InputFile> diffExpressionFiles = diffExpression.stream().collect(HashMap::new, (map, file) -> {
            String name = file.getName();
            map.put(name.substring(0, name.indexOf(".")), addInput(inputDiffExpression, file));
        }, HashMap::putAll);

        affinityRatios.forEach((pairing, hmMap) -> {
            InputFile diffExpressionFile = diffExpressionFiles.get(pairing);
            OutputFile inPairingAffinityRatios = new OutputFile(inputAffinityRatios, pairing);
            OutputFile outPairing = new OutputFile(outputDirectory, pairing);

            outputFiles.computeIfAbsent(pairing, s -> new HashMap<>());

            hmMap.forEach((hm, affinityRatioFile) -> {
                OutputFile outHm = addOutput(outPairing, hm + ".tsv");
                InputFile affinityRatioInput = addInput(inPairingAffinityRatios, affinityRatioFile);
                bridge.put(outHm, Pair.of(affinityRatioInput, diffExpressionFile));
                outputFiles.get(pairing).put(hm, outHm);
            });
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((outputFile, pair) -> add(() -> {
                List<String> commandElements =
                        List.of("python3", script.getAbsolutePath(), pair.getLeft().getAbsolutePath(),
                                pair.getRight().getAbsolutePath(), outputFile.getAbsolutePath(),
                                // Gene ID column index (inside diffExpression file)
                                "--geneIDs", "0",
                                // Log2FC column index (inside diffExpression file)
                                "--expressionC", "1");

                String command = String.join(" ", commandElements);

                executeAndWait(command, true);

                return true;
            }));
        }};
    }
}
