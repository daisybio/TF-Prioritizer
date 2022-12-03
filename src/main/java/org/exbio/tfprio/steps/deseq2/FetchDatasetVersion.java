package org.exbio.tfprio.steps.deseq2;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.Collection;
import java.util.HashSet;
import java.util.concurrent.Callable;

public class FetchDatasetVersion extends ExecutableStep {
    public final OutputFile outputFile = addOutput("version.tsv");
    private final RequiredConfig<String> species = new RequiredConfig<>(Configs.deSeq2.species);

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                URL url = new URL("http://www.ensembl.org/biomart/martservice?type=datasets&mart=ENSEMBL_MART_ENSEMBL");
                HttpURLConnection conn = (HttpURLConnection) url.openConnection();
                conn.setRequestMethod("GET");
                try (BufferedReader reader = new BufferedReader(new InputStreamReader(conn.getInputStream()));
                     BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
                    for (String line; (line = reader.readLine()) != null; ) {
                        if (line.contains(species.get())) {
                            String[] split = line.split("\t");
                            writer.write(String.join("\t", split[1], split[2], split[4]));
                            return true;
                        }
                    }
                }
                logger.error("Could not find dataset version for species " + species.get());
                return false;
            });
        }};
    }
}
