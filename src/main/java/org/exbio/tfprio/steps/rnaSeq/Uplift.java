package org.exbio.tfprio.steps.rnaSeq;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.util.Collection;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.FileManagement.copyFile;
import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;
import static org.exbio.tfprio.util.Helpers.getDatasetVersion;

public class Uplift extends ExecutableStep {
    public final OutputFile outputFile;
    private final RequiredConfig<String> species = new RequiredConfig<>(Configs.deSeq2.speciesBiomart);
    private final RequiredConfig<String> speciesRefGenome = new RequiredConfig<>(Configs.deSeq2.speciesRefGenome);
    private final InputFile genePositions;

    private final InputFile script;

    private final RequiredConfig<Map<String, String>> grcSynonymDict =
            new RequiredConfig<>(Configs.internalConfigs.grcSynonymDict);

    public Uplift(OutputFile genePositions) {
        super(false, genePositions);

        this.genePositions = addInput(genePositions);
        outputFile = addOutput("uplifted.tsv");

        script = addInput(getClass().getResourceAsStream("uplift.py"), "uplift.py");
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                String biomartVersion = getDatasetVersion(species.get());

                if (biomartVersion == null) {
                    logger.error("Could not find biomart version for species " + species.get());
                    return false;
                }

                if (!grcSynonymDict.get().containsKey(biomartVersion)) {
                    logger.error("Could not find synonym for biomart version " + biomartVersion);
                    return false;
                }

                String originalVersion = grcSynonymDict.get().get(biomartVersion);
                String refGenomeVersion = speciesRefGenome.get();

                if (originalVersion.equals(refGenomeVersion)) {
                    logger.info("No need to uplift, versions match.");
                    copyFile(genePositions, outputFile);
                    return true;
                }

                executeAndWait("python3 " + script + " " +
                        String.join(" ", genePositions.getAbsolutePath(), outputFile.getAbsolutePath(), originalVersion,
                                refGenomeVersion, "2", "3", "4"), true);
                return true;
            });
        }};
    }
}
