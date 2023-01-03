package org.exbio.tfprio.steps.logos;

import org.apache.commons.lang3.tuple.Pair;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;

import java.io.File;
import java.util.*;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class PredictedLogo extends ExecutableStep {
    public final Map<String, OutputFile> outputFiles = new HashMap<>();
    private final Map<InputFile, OutputFile> bridge = new HashMap<>();
    private final InputFile script;

    public PredictedLogo(Map<String, OutputFile> hmDirectories) {
        super(false, hmDirectories.values());

        script = addInput(getClass().getResourceAsStream("predicted.R"), "predicted.R");

        hmDirectories.forEach((hm, directory) -> {
            InputFile inputHm = addInput(directory);
            OutputFile outputHm = addOutput(hm);
            outputFiles.put(hm, outputHm);
            bridge.put(inputHm, outputHm);
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((inputHm, outputHm) -> {
                File[] tfFiles = inputHm.listFiles();

                if (tfFiles == null) {
                    logger.warn("No files found in directory: " + inputHm);
                    return;
                }

                Arrays.stream(tfFiles).map(tfFile -> {
                    String name = tfFile.getName();
                    return Pair.of(name.substring(0, name.lastIndexOf(".")), tfFile);
                }).forEach(tfPair -> {
                    add(() -> {
                        String tfName = tfPair.getKey();

                        File tfFile = tfPair.getValue();
                        File logoFile = new File(outputHm, tfName + ".png");
                        File outputFile = new File(outputHm, tfName + ".tsv");

                        String command =
                                String.format("Rscript %s --input_file %s --plot_file %s --output_file %s ", script,
                                        tfFile, logoFile, outputFile);

                        executeAndWait(command, true);

                        return true;
                    });
                });
            });
        }};
    }
}
