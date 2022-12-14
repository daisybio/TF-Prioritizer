package org.exbio.tfprio.steps.TEPIC;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;

import java.io.File;
import java.net.URL;
import java.util.Collection;
import java.util.HashSet;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.FileManagement.*;

public class FetchTEPIC extends ExecutableStep {
    public final OutputFile outputFile = addOutput("TEPIC");

    public FetchTEPIC() {
        super();
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                deleteFileStructure(outputFile);
                String path = "/org/exbio/tfprio/steps/TEPIC";
                URL resource = getClass().getResource(path);
                if (resource == null || !new File(resource.getFile()).exists()) {
                    // Usually entered when running via Jar
                    logger.info("Jar mode");
                    copyResources(path, outputFile.toPath());
                } else {
                    // Usually entered when running via IDE
                    logger.info("IDE mode");
                    copyDirectory(new File(resource.getFile()), outputFile, false);
                }
                makeAllChildrenExecutable(outputFile);
                return true;
            });
        }};
    }

    @Override
    protected boolean mayBeSkipped() {
        return false;
    }
}
