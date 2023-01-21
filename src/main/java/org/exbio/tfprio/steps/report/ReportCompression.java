package org.exbio.tfprio.steps.report;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;

import java.io.File;
import java.util.Collection;
import java.util.HashSet;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.FileManagement.writeFile;
import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class ReportCompression extends ExecutableStep {
    public final OutputFile compressed = addOutput("report.zip");
    private final InputFile uncompressed;

    public ReportCompression(OutputFile reportDirectory) {
        super(false, reportDirectory);

        uncompressed = addInput(reportDirectory);
    }

    @Override
    protected boolean doCreateFiles() {
        return false;
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                String command =
                        "pushd " + uncompressed.getAbsolutePath() + " && zip -qr " + compressed.getAbsolutePath() +
                                " ./* ";
                File scriptFile = new File(outputDirectory, "compress.sh");
                writeFile(scriptFile, command);

                executeAndWait("bash " + scriptFile.getAbsolutePath(), true);
                
                return true;
            });
        }};
    }
}
