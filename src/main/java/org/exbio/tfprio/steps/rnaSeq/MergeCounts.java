package org.exbio.tfprio.steps.rnaSeq;

import org.apache.commons.lang3.tuple.Pair;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.pipejar.util.FileFilters.Filters;
import org.exbio.tfprio.configs.Configs;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.exbio.pipejar.util.FileManagement.readLines;

public class MergeCounts extends ExecutableStep<Configs> {
    public final Map<String, OutputFile> outputFiles = new HashMap<>();
    private final RequiredConfig<File> rnaSeqDirectory = new RequiredConfig<>(configs.inputConfigs.rnaSeq);
    private final InputFile geneIDFile;
    private final Map<OutputFile, Collection<File>> bridge = new HashMap<>();

    public MergeCounts(Configs configs, OutputFile ensgFile) {
        super(configs, false, ensgFile);
        geneIDFile = addInput(ensgFile);
        // Register all input files and their corresponding output files
        InputFile input = addInput(rnaSeqDirectory);
        Arrays.stream(Objects.requireNonNull(input.listFiles(Filters.directoryFilter))).forEach(d_group -> {
            OutputFile f_groupOut = addOutput(d_group.getName() + ".tsv");
            outputFiles.put(d_group.getName(), f_groupOut);
            bridge.put(f_groupOut, new ArrayList<>() {{
                this.addAll(Arrays.asList(Objects.requireNonNull(d_group.listFiles(Filters.fileFilter))));
            }});
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            try {
                // Read gene IDs
                List<String> geneIDs = readLines(geneIDFile);

                bridge.forEach((outputFile, inputFiles) -> add(() -> {
                    // Read counts from all input files
                    TreeMap<String, List<Integer>> fileCounts =
                            inputFiles.stream().collect(TreeMap::new, (map, file) -> {
                                try {
                                    List<String> content = readLines(file);
                                    String fileName = file.getName();
                                    String colName = fileName.substring(0, fileName.lastIndexOf('.'));
                                    colName = colName.replace("-", "_");
                                    map.put(colName, content.stream().map(Double::parseDouble).map(
                                            d -> (int) Math.round(d)).collect(Collectors.toList()));
                                } catch (IOException e) {
                                    e.printStackTrace();
                                }
                            }, TreeMap::putAll);

                    // Check if sizes match
                    if (!fileCounts.values().stream().allMatch(l -> l.size() == geneIDs.size())) {
                        throw new RuntimeException("The number of gene IDs is " + geneIDs.size() +
                                ". The sizes of the following files differ: " +
                                fileCounts.entrySet().stream().filter(e -> e.getValue().size() != geneIDs.size()).map(
                                        Map.Entry::getKey).collect(Collectors.joining(", ")));
                    }

                    Map<String, List<Integer>> geneIdCounts = IntStream.range(0, geneIDs.size()).mapToObj(i -> {
                        String geneID = geneIDs.get(i);
                        return Pair.of(geneID, fileCounts.navigableKeySet().stream().map(
                                colName -> fileCounts.get(colName).get(i)).collect(Collectors.toList()));
                    }).collect(Collectors.toMap(Pair::getKey, Pair::getValue,
                            (old, next) -> IntStream.range(0, old.size()).mapToObj(
                                    i -> old.get(i) + next.get(i)).collect(Collectors.toList())));

                    // Write output file
                    try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
                        // Write header
                        writer.write("gene_id\t" + String.join("\t", fileCounts.navigableKeySet()));
                        writer.newLine();
                        // For all the genes
                        geneIdCounts.entrySet().stream().filter(geneIdEntry -> !geneIdEntry.getKey().equals("-"))
                                    // Map gene index to an output line
                                    .map(i -> i.getKey() + "\t" + String.join("\t", i.getValue().stream().map(
                                            Object::toString).toList())).sorted().forEachOrdered(line -> {
                                        try {
                                            writer.write(line);
                                            writer.newLine();
                                        } catch (IOException e) {
                                            throw new RuntimeException(e);
                                        }
                                    });
                    }
                    return true;
                }));
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }};
    }
}
