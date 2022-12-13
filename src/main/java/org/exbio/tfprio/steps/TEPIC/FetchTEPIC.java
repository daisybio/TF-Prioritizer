package org.exbio.tfprio.steps.TEPIC;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;

import java.io.File;
import java.net.URL;
import java.util.Collection;
import java.util.HashSet;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.FileManagement.copyDirectory;
import static org.exbio.pipejar.util.FileManagement.makeAllChildrenExecutable;

public class FetchTEPIC extends ExecutableStep {
    public final OutputFile outputFile = addOutput("TEPIC");

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                URL resource = (getClass().getResource("/org/exbio/tfprio/steps/TEPIC"));
                if (resource == null) {
                    copyResources("/org/exbio/tfprio/steps/TEPIC", outputFile.toPath());
                } else {
                    copyDirectory(new File(resource.getFile()), outputFile, false);
                }
                makeAllChildrenExecutable(outputFile);
                return true;
            });
        }};
    }
}
