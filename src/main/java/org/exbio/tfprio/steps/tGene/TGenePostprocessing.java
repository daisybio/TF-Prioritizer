package org.exbio.tfprio.steps.tGene;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.OptionalConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;
import org.exbio.tfprio.lib.GeneRegion;
import org.exbio.tfprio.lib.MultiGeneRegion;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import static org.exbio.pipejar.util.FileManagement.readLines;

public class TGenePostprocessing extends ExecutableStep {
    public final Map<String, Map<String, OutputFile>> outputFiles = new HashMap<>();
    private final Map<String, InputFile> tpmFiles = new HashMap<>();
    private final Map<String, Map<String, Collection<InputFile>>> linkFiles = new HashMap<>();
    private final OptionalConfig<Double> tpmFilter = new OptionalConfig<>(Configs.deSeq2.tpmFilter, false);

    public TGenePostprocessing(Map<String, OutputFile> tpmFiles,
                               Map<String, Map<String, Collection<OutputFile>>> linkFiles) {
        super(false, tpmFiles.values(),
                linkFiles.values().stream().flatMap(map -> map.values().stream()).flatMap(Collection::stream).toArray(
                        OutputFile[]::new));

        tpmFiles.forEach((group, file) -> this.tpmFiles.put(group, addInput(file)));

        linkFiles.forEach((group, hmMap) -> {
            this.linkFiles.put(group, new HashMap<>());
            outputFiles.put(group, new HashMap<>());
            hmMap.forEach((hm, sampleFiles) -> {
                this.linkFiles.get(group).put(hm, new HashSet<>());
                outputFiles.get(group).put(hm, addOutput("links_" + group + "_" + hm + ".tsv"));
                sampleFiles.forEach(sampleFile -> this.linkFiles.get(group).get(hm).add(addInput(sampleFile)));
            });
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            linkFiles.forEach((group, hmMap) -> {
                final TreeSet<String> validGenes;
                try {
                    validGenes = readLines(tpmFiles.get(group)).stream().skip(1).map(line -> line.split("\t")).filter(
                            split -> !tpmFilter.isSet() || Double.parseDouble(split[1]) > tpmFilter.get()).map(
                            split -> split[0]).collect(Collectors.toCollection(TreeSet::new));
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }

                hmMap.forEach((hm, sampleFiles) -> add(() -> {
                    Map<InputFile, List<GeneRegion>> sampleRegions = new HashMap<>();

                    sampleFiles.forEach(sampleFile -> {
                        try {
                            List<GeneRegion> currentRegions = GeneRegion.removeGeneRegionDuplicates(
                                    readLines(sampleFile).stream().skip(1).map(line -> line.split("\t")).filter(
                                            split -> split.length > 6).map(split -> {
                                        String gene = split[0].substring(0, split[0].indexOf('.'));
                                        String location = split[6];

                                        String chromosome = location.substring(0, location.indexOf(":"));
                                        int start = Integer.parseInt(
                                                location.substring(location.indexOf(":") + 1, location.indexOf("-")));
                                        int end = Integer.parseInt(location.substring(location.indexOf("-") + 1));

                                        return new GeneRegion(chromosome, start, end, gene);
                                    }).filter(region -> validGenes.contains(region.getId())).sorted().toList(), true);

                            sampleRegions.put(sampleFile, currentRegions);
                        } catch (IOException e) {
                            throw new RuntimeException(e);
                        }
                    });


                    List<MultiGeneRegion> combined = GeneRegion.mergeSameRegions(GeneRegion.removeGeneRegionDuplicates(
                            sampleRegions.values().stream().flatMap(Collection::stream).sorted().toList(), true));

                    try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFiles.get(group).get(hm)))) {
                        writer.write("CHR\tLEFT_BORDER\tRIGHT_BORDER\tENSGS");
                        writer.newLine();

                        combined.forEach(region -> {
                            try {
                                writer.write(region.toString());
                                writer.newLine();
                            } catch (IOException e) {
                                throw new RuntimeException(e);
                            }
                        });
                    }
                    return true;
                }));
            });
        }};
    }
}
