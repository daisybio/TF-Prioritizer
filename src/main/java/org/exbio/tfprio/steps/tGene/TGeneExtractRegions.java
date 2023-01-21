package org.exbio.tfprio.steps.tGene;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;
import org.exbio.tfprio.lib.GeneRegion;

import java.io.*;
import java.util.Collection;
import java.util.HashSet;
import java.util.TreeSet;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

public class TGeneExtractRegions extends ExecutableStep<Configs> {
    public final OutputFile outputFile;
    private final InputFile inputFile;

    public TGeneExtractRegions(Configs configs, OutputFile input) {
        super(configs, false, input);
        this.inputFile = addInput(input);
        this.outputFile = addOutput("regions.tsv");
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                TreeSet<GeneRegion> regions = new TreeSet<>();

                try (BufferedReader reader = new BufferedReader(new FileReader(inputFile))) {
                    for (String line = reader.readLine(); line != null; line = reader.readLine()) {
                        String[] fields = line.split("\t");
                        String chr = fields[0];
                        int start = Integer.parseInt(fields[3]);
                        int end = Integer.parseInt(fields[4]);
                        String geneId = fields[8].split(";")[0].split(" ")[1].replace("\"", "");
                        geneId = geneId.substring(0, geneId.lastIndexOf("."));

                        regions.add(new GeneRegion(chr, start, end, geneId));
                    }
                }

                try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
                    regions.stream().collect(Collectors.groupingBy(GeneRegion::getId)).values().stream().map(
                            matchingRegions -> {
                                GeneRegion first = matchingRegions.get(0);
                                matchingRegions.stream().skip(1).forEach(first::merge);
                                return first;
                            }).sorted().forEachOrdered(mergedRegion -> {
                        try {
                            writer.write(mergedRegion.toString());
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
