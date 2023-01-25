package org.exbio.tfprio.steps.rnaSeq;

import org.apache.commons.lang3.tuple.Pair;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class SplitCounts extends ExecutableStep<Configs> {
    public final Map<String, OutputFile> outputFiles = new HashMap<>();
    private final InputFile geneIDFile;
    private final InputFile rnaSeqFile;
    private final RequiredConfig<File> rnaSeqConfig = new RequiredConfig<>(configs.inputConfigs.rnaSeq);
    private final Map<String, List<Integer>> stageColumns;

    public SplitCounts(Configs configs, OutputFile ensgFile) {
        super(configs, false, ensgFile);

        geneIDFile = addInput(ensgFile);
        rnaSeqFile = addInput(rnaSeqConfig);

        try (var reader = new BufferedReader(new FileReader(rnaSeqConfig.get()))) {
            String[] header = reader.readLine().split("\t");
            stageColumns = IntStream.range(2, header.length).mapToObj(i -> {
                String s = header[i];
                String stage = s.split("_")[0];

                return Pair.of(stage, i);
            }).collect(Collectors.groupingBy(Pair::getKey, Collectors.mapping(Pair::getValue, Collectors.toList())));
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }

        stageColumns.keySet().forEach(stage -> outputFiles.put(stage, addOutput(stage + ".tsv")));
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                try (var reader = new BufferedReader(new FileReader(rnaSeqFile))) {
                    Map<String, BufferedWriter> writers =
                            stageColumns.keySet().stream().collect(Collectors.toMap(stage -> stage, stage -> {
                                try {
                                    return new BufferedWriter(new FileWriter(outputFiles.get(stage)));
                                } catch (IOException e) {
                                    throw new UncheckedIOException(e);
                                }
                            }));

                    reader.lines().forEach(line -> {
                        String[] split = line.split("\t");
                        String geneID = split[0];

                        stageColumns.forEach((stage, columns) -> {
                            try {
                                writers.get(stage).write(geneID);
                                columns.forEach(column -> {
                                    try {
                                        writers.get(stage).write("\t" + split[column]);
                                    } catch (IOException e) {
                                        throw new UncheckedIOException(e);
                                    }
                                });
                                writers.get(stage).write("\n");
                            } catch (IOException e) {
                                throw new UncheckedIOException(e);
                            }
                        });
                    });

                    writers.values().forEach(writer -> {
                        try {
                            writer.close();
                        } catch (IOException e) {
                            throw new UncheckedIOException(e);
                        }
                    });
                }

                return true;
            });
        }};
    }
}
