package org.exbio.tfprio.steps.report;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;

import java.io.*;
import java.nio.file.Files;
import java.util.Collection;
import java.util.HashSet;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.FileManagement.writeFile;
import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class ReportCreation extends ExecutableStep {
    public final OutputFile outputFile = addOutput("report");
    private final InputFile reportDirectory;

    public ReportCreation(OutputFile reportDirectory) {
        super(false, reportDirectory);

        this.reportDirectory = addInput(reportDirectory);
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                String command = "pushd " + reportDirectory + " && npm install && ng build --output-path " +
                        outputFile.getAbsolutePath() + " && popd";

                File scriptFile = new File(outputFile, "create.sh");

                writeFile(scriptFile, command);

                executeAndWait("bash " + scriptFile.getAbsolutePath(), true);

                File index = new File(outputFile, "index.html");
                File processed = new File(outputFile, "processed_index.html");

                try (var reader = new BufferedReader(new FileReader(index));
                     var writer = new BufferedWriter(new FileWriter(processed))) {
                    reader.lines().map(line -> line.replace("type=\"module\"", "")).forEach(line -> {
                        try {
                            writer.write(line);
                            writer.newLine();
                        } catch (IOException e) {
                            throw new UncheckedIOException(e);
                        }
                    });
                }

                Files.delete(index.toPath());
                Files.move(processed.toPath(), index.toPath());

                return true;
            });
        }};
    }
}
