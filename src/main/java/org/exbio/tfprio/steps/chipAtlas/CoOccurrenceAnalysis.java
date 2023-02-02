package org.exbio.tfprio.steps.chipAtlas;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.util.Collection;
import java.util.HashSet;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class CoOccurrenceAnalysis extends ExecutableStep<Configs> {
    public final OutputFile outputFile = addOutput("coOccurrence.tsv");
    public final OutputFile merged = addOutput("merged.bed");
    private final InputFile inputDirectory;
    private final InputFile mergeScript;
    private final InputFile coOccurrenceScript;

    public CoOccurrenceAnalysis(Configs configs, OutputFile chipAtlasDirectory) {
        super(configs, false, chipAtlasDirectory);

        inputDirectory = addInput(chipAtlasDirectory);
        mergeScript = addInput(getClass().getResourceAsStream("merge.sh"), "merge.sh");
        coOccurrenceScript = addInput(getClass().getResourceAsStream("cooccurrence.py"), "cooccurrence.py");
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                String mergeCommand = "bash " + mergeScript + " " + inputDirectory + " " + merged;
                executeAndWait(mergeCommand, true);

                String coOccurrenceCommand = "python3 " + coOccurrenceScript + " " + merged + " " + outputFile;
                executeAndWait(coOccurrenceCommand, true);

                return true;
            });
        }};
    }
}
