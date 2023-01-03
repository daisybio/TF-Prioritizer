package org.exbio.tfprio.steps.logos;

import org.apache.commons.lang3.tuple.Pair;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.OptionalConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import static org.exbio.pipejar.util.FileManagement.readLines;

public class ExtractBindingSites extends ExecutableStep {
    public final Map<String, OutputFile> outputFiles = new HashMap<>();
    private final InputFile dcgFile;
    private final Map<OutputFile, Collection<InputFile>> bridge = new HashMap<>();
    private final OptionalConfig<Double> affinityCutoff = new OptionalConfig<>(Configs.tepic.affinityCutoff, false);

    public ExtractBindingSites(OutputFile dcgFile,
                               Map<String, Map<String, Collection<OutputFile>>> groupHmSampleSequenceFiles) {
        super(false, groupHmSampleSequenceFiles.values().stream().flatMap(sub -> sub.values().stream()).flatMap(
                Collection::stream).toList(), dcgFile);
        this.dcgFile = addInput(dcgFile);

        Set<String> hms = groupHmSampleSequenceFiles.values().stream().flatMap(sub -> sub.keySet().stream()).collect(
                Collectors.toSet());

        hms.forEach(hm -> {
            OutputFile outputHm = addOutput(hm);
            outputFiles.put(hm, outputHm);
        });

        groupHmSampleSequenceFiles.forEach((group, hmMap) -> {
            OutputFile inputGroup = new OutputFile(inputDirectory, group);

            hmMap.forEach((hm, sampleDirs) -> {
                OutputFile inputHm = new OutputFile(inputGroup, hm);
                OutputFile outputHm = outputFiles.get(hm);

                sampleDirs.forEach(sampleDir -> {
                    InputFile inputSample = addInput(inputHm, sampleDir);
                    bridge.computeIfAbsent(outputHm, s -> new HashSet<>()).add(inputSample);
                });
            });
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            Set<String> tfNames;
            try {
                tfNames = readLines(dcgFile).stream().skip(1).map(line -> line.split("\t")[0]).collect(
                        Collectors.toSet());
            } catch (IOException e) {
                throw new RuntimeException(e);
            }

            bridge.forEach((output, inputs) -> add(() -> {
                Map<String, BufferedWriter> tfWriter = new HashMap<>() {{
                    tfNames.forEach(tf -> {
                        File tfFile = new File(output, tf + ".fa");
                        try {
                            put(tf, new BufferedWriter(new FileWriter(tfFile)));
                        } catch (IOException e) {
                            throw new RuntimeException(e);
                        }
                    });
                }};

                inputs.forEach(inputFile -> {
                    try (BufferedReader reader = new BufferedReader(new FileReader(inputFile))) {
                        reader.lines().filter(line -> !line.startsWith("#") && !line.startsWith("TF\tAFFINITY")).map(
                                line -> line.split("\t")).filter(split -> !affinityCutoff.isSet() ||
                                Double.parseDouble(split[1]) > affinityCutoff.get()).map(
                                split -> Pair.of(split[0], split[5])).filter(
                                pair -> tfWriter.containsKey(pair.getKey())).forEach(pair -> {
                            var writer = tfWriter.get(pair.getKey());
                            try {
                                writer.write(pair.getValue());
                                writer.newLine();
                            } catch (IOException e) {
                                throw new RuntimeException(e);
                            }
                        });
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                });

                tfWriter.values().forEach(writer -> {
                    try {
                        writer.close();
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                });


                return true;
            }));
        }};
    }
}
