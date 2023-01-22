package org.exbio.tfprio.steps.Dynamite;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.OptionalConfig;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.UsageConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class RunDynamite extends ExecutableStep<Configs> {
    public final Map<String, Map<String, OutputFile>> outputFiles = new HashMap<>();
    private final InputFile script;
    private final Map<InputFile, OutputFile> bridge = new HashMap<>();

    private final RequiredConfig<String> outVar = new RequiredConfig<>(configs.dynamite.outVar);
    private final OptionalConfig<Integer> oFolds = new OptionalConfig<>(configs.dynamite.oFolds, false);
    private final OptionalConfig<Integer> iFolds = new OptionalConfig<>(configs.dynamite.iFolds, false);
    private final OptionalConfig<Boolean> performance = new OptionalConfig<>(configs.dynamite.performance, false);
    private final OptionalConfig<Double> alpha = new OptionalConfig<>(configs.dynamite.alpha, false);
    private final OptionalConfig<Boolean> randomize = new OptionalConfig<>(configs.dynamite.randomize, false);

    public RunDynamite(Configs configs, Map<String, Map<String, OutputFile>> preparedData) {
        super(configs, false, preparedData.values().stream().flatMap(
                stringCollectionMap -> stringCollectionMap.values().stream()).toList());
        script = addInput(getClass().getResourceAsStream("/org/exbio/tfprio/steps/TEPIC/DYNAMITE/Scripts/DYNAMITE.R"),
                "DYNAMITE.R");

        preparedData.forEach((pairing, hmMap) -> {
            OutputFile inputPairing = new OutputFile(inputDirectory, pairing);
            OutputFile outputPairing = new OutputFile(outputDirectory, pairing);

            hmMap.forEach((hm, inDir) -> {
                InputFile inputDir = addInput(inputPairing, inDir);
                OutputFile outputDir = addOutput(outputPairing, hm);

                bridge.put(inputDir, outputDir);
                outputFiles.computeIfAbsent(pairing, s -> new HashMap<>()).put(hm, outputDir);
            });
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            Map<String, UsageConfig<?>> parameters = new HashMap<>() {{
                put("out_var", outVar);
                put("Ofolds", oFolds);
                put("Ifolds", iFolds);
                put("performance", performance);
                put("alpha", alpha);
                put("randomise", randomize);
            }};

            String parametersString = parameters.entrySet().stream().filter(entry -> entry.getValue().isSet()).map(
                    entry -> "--" + entry.getKey() + "=" + entry.getValue().get()).collect(Collectors.joining(" "));

            bridge.forEach((inDir, outDir) -> add(() -> {
                List<String> commandElements =
                        List.of("Rscript", script.getAbsolutePath(), "--dataDir=" + inDir.getAbsolutePath(),
                                "--outDir=" + outDir.getAbsolutePath());

                String command = String.join(" ", commandElements) + " " + parametersString;

                logger.debug(command);

                executeAndWait(command, true);

                return true;
            }));
        }};
    }
}
