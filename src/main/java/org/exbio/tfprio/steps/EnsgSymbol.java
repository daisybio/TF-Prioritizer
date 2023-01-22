package org.exbio.tfprio.steps;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import static org.exbio.pipejar.util.FileManagement.readLines;
import static org.exbio.tfprio.lib.Biomart.query;

public class EnsgSymbol extends ExecutableStep<Configs> {
    public final OutputFile outputFile = addOutput("ensgSymbol.tsv");
    private final InputFile geneIdsFile;
    private final RequiredConfig<String> species = new RequiredConfig<>(configs.deSeq2.speciesBiomart);


    public EnsgSymbol(Configs configs, OutputFile ensgFile) {
        super(configs, false, ensgFile);
        geneIdsFile = addInput(ensgFile);
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                List<String[]> result = query(species.get(), readLines(geneIdsFile).stream().skip(1).toList(),
                        List.of("ensembl_gene_id", "external_gene_name"));
                try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
                    result.forEach(line -> {
                        try {
                            writer.write(Arrays.stream(line).skip(1).map(String::toUpperCase).collect(
                                    Collectors.joining("\t")));
                            writer.newLine();
                        } catch (IOException e) {
                            e.printStackTrace();
                            throw new RuntimeException(e);
                        }
                    });
                }
                return true;
            });
        }};
    }
}
