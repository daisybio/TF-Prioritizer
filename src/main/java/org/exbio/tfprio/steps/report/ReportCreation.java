package org.exbio.tfprio.steps.report;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.*;
import java.nio.file.Files;
import java.util.Collection;
import java.util.HashSet;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.FileManagement.copyDirectory;
import static org.exbio.pipejar.util.FileManagement.writeFile;
import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class ReportCreation extends ExecutableStep<Configs> {
    public final OutputFile outputFile = addOutput("report");
    private final InputFile reportDirectory;

    public ReportCreation(Configs configs, OutputFile reportDirectory) {
        super(configs, false, reportDirectory);

        this.reportDirectory = addInput(reportDirectory);
    }

    @Override
    protected boolean doCreateFiles() {
        return false;
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                // The temp directories are a workaround for the problems that occurred when running the
                // report generation on the server. The exact problem was that somehow the generation of the report
                // tried to make sure that the output directory is empty, but it failed to do so. A possible reason
                // for this are permissions, since everything else should be the same due to docker.
                // TODO: Find a better solution for this problem.
                File tempAngular = new File("/srv/temp/angular");
                File tempOutput = new File("/srv/temp/output");

                copyDirectory(reportDirectory, tempAngular);

                String command = "pushd " + tempAngular + " && npm install && ng build --output-path " +
                        tempOutput.getAbsolutePath() + " && popd";

                File scriptFile = new File(outputDirectory, "create.sh");

                writeFile(scriptFile, command);

                executeAndWait("bash " + scriptFile.getAbsolutePath(), true);

                File index = new File(tempOutput, "index.html");
                File processed = new File(tempOutput, "processed_index.html");

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

                copyDirectory(tempOutput, outputFile);

                return true;
            });
        }};
    }
}
