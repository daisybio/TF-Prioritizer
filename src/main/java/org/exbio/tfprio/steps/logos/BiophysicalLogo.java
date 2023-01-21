package org.exbio.tfprio.steps.logos;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.File;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.Objects;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class BiophysicalLogo extends ExecutableStep<Configs> {
    public final OutputFile outputFile = addOutput("logos");
    private final InputFile biophysicalModelDirectory;
    private final InputFile script;

    public BiophysicalLogo(Configs configs, OutputFile biophysicalModelDirectory) {
        super(configs, false, biophysicalModelDirectory);

        this.biophysicalModelDirectory = addInput(biophysicalModelDirectory);
        script = addInput(getClass().getResourceAsStream("biophysical.py"), "biophysical.py");
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            Arrays.stream(Objects.requireNonNull(biophysicalModelDirectory.listFiles(File::isFile))).forEach(
                    inputFile -> add(() -> {
                        String inputName = inputFile.getName();
                        File output = new File(outputFile, inputName.substring(0, inputName.lastIndexOf(".")) + ".png");

                        String command = String.format("python3 %s %s %s", script, inputFile, output);
                        executeAndWait(command, true);

                        return true;
                    }));
        }};
    }
}
