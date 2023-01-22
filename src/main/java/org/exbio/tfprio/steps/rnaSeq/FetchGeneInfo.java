package org.exbio.tfprio.steps.rnaSeq;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.FileManagement.readLines;

public class FetchGeneInfo extends ExecutableStep<Configs> {
    public final OutputFile outputFile = addOutput("geneInfo.tsv");
    private final InputFile geneIdsFile;
    private final RequiredConfig<String> species = new RequiredConfig<>(configs.deSeq2.speciesBiomart);

    public FetchGeneInfo(Configs configs, OutputFile ensgFile) {
        super(configs, false, ensgFile);

        geneIdsFile = addInput(ensgFile);
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                List<String> geneIds = readLines(geneIdsFile).stream().skip(1).filter(id -> !id.equals("-")).toList();

                org.exbio.tfprio.lib.FetchGeneInfo fetch =
                        new org.exbio.tfprio.lib.FetchGeneInfo(species.get(), geneIds, "ensembl_gene_id");

                try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
                    fetch.getLines().forEach(line -> {
                        try {
                            writer.write(line);
                            writer.newLine();
                        } catch (IOException e) {
                            throw new RuntimeException(e);
                        }
                    });
                }
                return true;
            });
        }};
    }
}
