package org.exbio.tfprio.steps.rnaSeq;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;
import org.exbio.tfprio.lib.Biomart;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.exbio.pipejar.util.FileManagement.readLines;

public class FetchGeneInfo extends ExecutableStep {
    private static final int batchSize = 1000;
    private final RequiredConfig<File> geneIdsFile = new RequiredConfig<>(Configs.inputConfigs.geneIDs);
    private final RequiredConfig<String> species = new RequiredConfig<>(Configs.deSeq2.speciesBiomart);
    private final Map<OutputFile, List<String>> bridge = new HashMap<>();

    public FetchGeneInfo() {
        super();

        InputFile inputFile = addInput(geneIdsFile);

        try {
            List<String> geneIds = readLines(inputFile).stream().skip(1).toList();

            IntStream.range(0, geneIds.size()).boxed().collect(Collectors.groupingBy(i -> i / batchSize)).forEach(
                    (batchId, indices) -> {
                        OutputFile f_out = addOutput("geneInfo_" + batchId + ".tsv");
                        bridge.put(f_out, indices.stream().map(geneIds::get).toList());
                    });
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((outputFile, geneIds) -> add(() -> {
                Biomart mart = new Biomart(species.get(), geneIds);
                try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
                    mart.getLines().forEach(line -> {
                        try {
                            writer.write(line);
                            writer.newLine();
                        } catch (IOException e) {
                            e.printStackTrace();
                            throw new RuntimeException(e);
                        }
                    });
                }
                return true;
            }));
        }};
    }
}
