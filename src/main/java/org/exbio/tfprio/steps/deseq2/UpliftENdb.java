package org.exbio.tfprio.steps.deseq2;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Collection;
import java.util.HashSet;
import java.util.Objects;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.FileManagement.softLink;
import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class UpliftENdb extends ExecutableStep {

    public final OutputFile outputFile;
    private final InputFile inputFile;
    private final String scriptPath;

    private final RequiredConfig<String> speciesRefGenome = new RequiredConfig<>(Configs.deSeq2.speciesRefGenome);

    public UpliftENdb(OutputFile enDBfile) {
        super(false, enDBfile);

        inputFile = addInput(enDBfile);
        outputFile = addOutput("upliftedENdb.tsv");
        scriptPath = Objects.requireNonNull(getClass().getResource("uplift.py")).getPath();
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                final String genome;
                try (BufferedReader reader = new BufferedReader(new FileReader(inputFile))) {
                    String line = reader.readLine();
                    genome = line.split("\t")[4];
                }

                if (genome.equals(speciesRefGenome.get())) {
                    logger.info("No need to uplift, versions match.");
                    softLink(outputFile, inputFile);
                    return true;
                }

                executeAndWait("python3 " + scriptPath + " " +
                        String.join(" ", inputFile.getPath(), outputFile.getPath(), genome, speciesRefGenome.get(), "0",
                                "1", "2"), true);

                return true;
            });
        }};
    }

    @Override
    protected boolean mayBeSkipped() {
        return false;
    }
}
