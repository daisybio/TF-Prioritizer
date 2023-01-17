package org.exbio.tfprio.steps.report;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;

import java.util.Collection;
import java.util.HashSet;
import java.util.concurrent.Callable;

public class ReportPostprocessing extends ExecutableStep {
    public final OutputFile outputFile = addOutput("report");
    private final InputFile inputDirectory;

    public ReportPostprocessing(OutputFile reportDirectory) {
        super(false, reportDirectory);

        this.inputDirectory = addInput(reportDirectory);
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                

                return true;
            });
        }};
    }
}
