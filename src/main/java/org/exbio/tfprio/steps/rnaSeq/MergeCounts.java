package org.exbio.tfprio.steps.rnaSeq;

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

public class MergeCounts extends ExecutableStep {
    public final Map<String, OutputFile> outputFiles = new HashMap<>();
    private final RequiredConfig<File> rnaSeqDirectory = new RequiredConfig<>(Configs.inputConfigs.rnaSeq);
    private final InputFile geneIDFile;
    private final Map<OutputFile, Collection<File>> bridge = new HashMap<>();

    public MergeCounts(OutputFile ensgFile) {
        super(false, ensgFile);
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
                    TreeMap<String, List<Integer>> counts = inputFiles.stream().collect(TreeMap::new, (map, file) -> {
                        try {
                            List<String> content = readLines(file);
                            String fileName = file.getName();
                            String colName = fileName.substring(0, fileName.lastIndexOf('.'));
                            colName = colName.replace("-", "_");
                            map.put(colName,
                                    content.stream().map(Double::parseDouble).map(d -> (int) Math.round(d)).collect(
                                            Collectors.toList()));
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                    }, TreeMap::putAll);

                    // Check if sizes match
                    if (!counts.values().stream().allMatch(l -> l.size() == geneIDs.size())) {
                        throw new RuntimeException("The number of gene IDS is " + geneIDs.size() +
                                ". The sizes of the following files differ: " +
                                counts.entrySet().stream().filter(e -> e.getValue().size() != geneIDs.size()).map(
                                        Map.Entry::getKey).collect(Collectors.joining(", ")));
                    }

                    // Write output file
                    try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
                        // Write header
                        writer.write("gene_id\t" + String.join("\t", counts.navigableKeySet()));
                        writer.newLine();
                        // For all the genes
                        IntStream.range(0, geneIDs.size()).filter(i -> !geneIDs.get(i).equals("-"))
                                 // Map gene index to an output line
                                 .mapToObj(i -> {
                                     List<Integer> countsAtI = counts.values().stream().map(l -> l.get(i)).toList();
                                     return geneIDs.get(i) + "\t" +
                                             String.join("\t", countsAtI.stream().map(Object::toString).toList());
                                 }).sorted().forEachOrdered(line -> {
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
