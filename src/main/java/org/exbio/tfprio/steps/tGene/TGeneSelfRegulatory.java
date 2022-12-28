package org.exbio.tfprio.steps.tGene;

import org.exbio.pipejar.pipeline.ExecutableStep;

import java.util.Collection;
import java.util.HashSet;
import java.util.concurrent.Callable;

public class TGeneSelfRegulatory extends ExecutableStep {
    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                logger.error("Not implemented yet");
                return false;
            });
        }};
    }
}
