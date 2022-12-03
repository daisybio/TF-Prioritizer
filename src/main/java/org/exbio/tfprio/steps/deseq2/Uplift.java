package org.exbio.tfprio.steps.deseq2;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.Collection;
import java.util.HashSet;
import java.util.Map;
import java.util.Objects;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.FileManagement.copyFile;
import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class Uplift extends ExecutableStep {
    public final OutputFile outputFile;
    private final RequiredConfig<String> species = new RequiredConfig<>(Configs.deSeq2.species);
    private final RequiredConfig<String> speciesRefGenome = new RequiredConfig<>(Configs.deSeq2.speciesRefGenome);
    private final InputFile genePositions;

    private final String scriptPath;

    private final RequiredConfig<Map<String, String>> grcSynonymDict =
            new RequiredConfig<>(Configs.internalConfigs.grcSynonymDict);

    public Uplift(OutputFile genePositions) {
        super(false, genePositions);

        this.genePositions = addInput(genePositions);
        outputFile = addOutput("uplifted.tsv");

        scriptPath = Objects.requireNonNull(getClass().getResource("uplift.py")).getPath();
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

                executeAndWait("python3 " + scriptPath + " " +
                        String.join(" ", genePositions.getAbsolutePath(), outputFile.getAbsolutePath(), originalVersion,
                                refGenomeVersion), true);
                return true;
            });
        }};
    }

    @Override
    protected boolean mayBeSkipped() {
        return false;
    }

    private String getDatasetVersion(String species) throws IOException {
        URL url = new URL("http://www.ensembl.org/biomart/martservice?type=datasets&mart=ENSEMBL_MART_ENSEMBL");
        HttpURLConnection conn = (HttpURLConnection) url.openConnection();
        conn.setRequestMethod("GET");
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(conn.getInputStream()))) {
            for (String line; (line = reader.readLine()) != null; ) {
                if (line.contains(species)) {
                    String[] split = line.split("\t");
                    return split[4];
                }
            }
        }
        return null;
    }
}
