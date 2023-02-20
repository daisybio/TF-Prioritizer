package org.exbio.tfprio.steps.report;

import org.apache.commons.lang3.tuple.Pair;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;
import org.json.JSONArray;
import org.json.JSONObject;

import java.io.*;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static java.util.stream.Collectors.*;
import static org.exbio.pipejar.util.FileManagement.*;

public class ReportPreprocessing extends ExecutableStep<Configs> {
    public final OutputFile outputFile = addOutput("report");
    private final InputFile ensgSymbol;
    private final InputFile angularInput;
    private final Map<String, InputFile> hmDcgFiles = new HashMap<>();
    private final Map<String, InputFile> pairingDeseq = new HashMap<>();
    private final Map<String, InputFile> groupMeanExpression = new HashMap<>();
    private final Map<String, InputFile> groupTpm = new HashMap<>();
    private final Map<String, Map<String, InputFile>> pairingHmHeatmapDir = new HashMap<>();
    private final Map<String, Map<String, InputFile>> pairingHmIgvDir = new HashMap<>();
    private final Map<String, InputFile> hmDistributionPlots = new HashMap<>();
    private final Map<String, Map<String, InputFile>> pairingHmRegressionCoefficients = new HashMap<>();
    private final InputFile biophysicalLogo;
    private final InputFile tfSequence;
    private final InputFile coOccurrence;
    private final Map<String, InputFile> importantLoci = new HashMap<>();
    private final Map<String, Pair<InputFile, InputFile>> topLog2fc = new HashMap<>();

    public ReportPreprocessing(Configs configs, OutputFile ensgSymbol, Map<String, OutputFile> hmDcgFiles,
                               Map<String, OutputFile> deseqResults, Map<String, OutputFile> groupMeanExpression,
                               Map<String, OutputFile> groupTpm, OutputFile biophysicalLogo, OutputFile tfSequence,
                               Map<String, Map<String, OutputFile>> pairingHmHeatmapDir,
                               Map<String, Map<String, OutputFile>> pairingHmIgvDir,
                               Map<String, OutputFile> hmDistributionPlots,
                               Map<String, Map<String, OutputFile>> pairingHmRegressionCoefficients,
                               OutputFile coOccurrence, Map<String, OutputFile> importantLoci,
                               Map<String, Pair<OutputFile, OutputFile>> topLog2fc) {
        super(configs, false, hmDcgFiles.values(), deseqResults.values(), groupMeanExpression.values(),
                groupTpm.values(), List.of(ensgSymbol), hmDistributionPlots.values(),
                pairingHmHeatmapDir.values().stream().flatMap(m -> m.values().stream()).collect(Collectors.toList()),
                pairingHmIgvDir.values().stream().flatMap(m -> m.values().stream()).collect(Collectors.toList()),
                pairingHmRegressionCoefficients.values().stream().flatMap(m -> m.values().stream()).collect(
                        Collectors.toList()));

        this.ensgSymbol = addInput(ensgSymbol);
        this.biophysicalLogo = addInput(biophysicalLogo);
        this.tfSequence = addInput(tfSequence);
        if (coOccurrence != null) {
            this.coOccurrence = addInput(coOccurrence);
        } else {
            this.coOccurrence = null;
        }
        OutputFile rankDir = new OutputFile(inputDirectory, "ranks");
        hmDcgFiles.forEach((hm, rankFile) -> this.hmDcgFiles.put(hm, addInput(rankDir, rankFile)));

        OutputFile deseqDir = new OutputFile(inputDirectory, "deseq");
        deseqResults.forEach((pairing, deseqFile) -> this.pairingDeseq.put(pairing, addInput(deseqDir, deseqFile)));

        OutputFile meanDir = new OutputFile(inputDirectory, "meanExpression");
        groupMeanExpression.forEach(
                (group, meanFile) -> this.groupMeanExpression.put(group, addInput(meanDir, meanFile)));

        OutputFile tpmDir = new OutputFile(inputDirectory, "tpm");
        groupTpm.forEach((group, tpmFile) -> this.groupTpm.put(group, addInput(tpmDir, tpmFile)));

        OutputFile heatmapDir = new OutputFile(inputDirectory, "heatmaps");
        pairingHmHeatmapDir.forEach((pairing, hmDirectories) -> {
            OutputFile pairingDir = new OutputFile(heatmapDir, pairing);
            hmDirectories.forEach((hm, hmDir) -> {
                InputFile input = addInput(pairingDir, hmDir);
                this.pairingHmHeatmapDir.computeIfAbsent(pairing, k -> new HashMap<>()).put(hm, input);
            });
        });

        OutputFile importantLociDir = new OutputFile(inputDirectory, "importantLoci");
        importantLoci.forEach((pairing, importantLociFile) -> this.importantLoci.put(pairing,
                addInput(importantLociDir, importantLociFile)));

        OutputFile topLog2fcDir = new OutputFile(inputDirectory, "topLog2fc");
        topLog2fc.forEach((pairing, topLog2fcFiles) -> {
            OutputFile pairingDir = new OutputFile(topLog2fcDir, pairing);
            this.topLog2fc.put(pairing, Pair.of(addInput(pairingDir, topLog2fcFiles.getLeft()),
                    addInput(pairingDir, topLog2fcFiles.getRight())));
        });

        OutputFile igvDir = new OutputFile(inputDirectory, "igv");
        pairingHmIgvDir.forEach((pairing, hmDirectories) -> {
            OutputFile pairingDir = new OutputFile(igvDir, pairing);
            hmDirectories.forEach((hm, hmDir) -> {
                InputFile input = addInput(pairingDir, hmDir);
                this.pairingHmIgvDir.computeIfAbsent(pairing, k -> new HashMap<>()).put(hm, input);
            });
        });

        OutputFile distributionDir = new OutputFile(inputDirectory, "distributionPlots");
        hmDistributionPlots.forEach((hm, distributionPlot) -> this.hmDistributionPlots.put(hm,
                addInput(distributionDir, distributionPlot)));

        OutputFile thresholdDir = new OutputFile(inputDirectory, "thresholdPlots");
        pairingHmRegressionCoefficients.forEach((pairing, hmFiles) -> {
            OutputFile pairingDir = new OutputFile(thresholdDir, pairing);
            hmFiles.forEach((hm, file) -> {
                InputFile input = addInput(pairingDir, file);
                this.pairingHmRegressionCoefficients.computeIfAbsent(pairing, k -> new HashMap<>()).put(hm, input);
            });
        });

        angularInput = new InputFile(inputDirectory, "report");

        String path = "/org/exbio/tfprio/steps/report/angular";
        URL resource = getClass().getResource(path);
        if (resource == null || !new File(resource.getFile()).exists()) {
            // Usually entered when running via Jar
            logger.info("Jar mode");
            try {
                copyResources(path, angularInput.toPath());
            } catch (URISyntaxException | IOException e) {
                throw new RuntimeException(e);
            }
        } else {
            // Usually entered when running via IDE
            logger.info("IDE mode");
            try {
                copyDirectory(new File(resource.getFile()), angularInput,
                        pathname -> !pathname.getName().endsWith(".class"));
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }

        logger.trace("Finished copying angular files");
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                copyDirectory(angularInput, outputFile);

                File srcDir = new File(outputFile, "src");

                Map<String, List<String>> symbolEnsgMap =
                        readLines(ensgSymbol).stream().map(line -> line.split("\t")).filter(
                                split -> split.length > 1).collect(
                                groupingBy(split -> split[1], mapping(split -> split[0], Collectors.toList())));

                logger.trace("Fetching tf symbols");

                Set<String> tfSymbols = hmDcgFiles.values().stream().flatMap(rankFile -> {
                    try {
                        return readLines(rankFile).stream().skip(1);
                    } catch (IOException e) {
                        throw new UncheckedIOException(e);
                    }
                }).map(line -> line.split("\t")[0]).collect(Collectors.toSet());

                logger.trace("Creating tf groups");

                Collection<TfGroup> groups = tfSymbols.stream().map(tfSymbol -> {
                    String[] tfParts = tfSymbol.split("::");

                    Collection<TranscriptionFactor> transcriptionFactors = Arrays.stream(tfParts).map(tfPart -> {
                        List<String> ensgs = symbolEnsgMap.get(tfPart);
                        return new TranscriptionFactor(tfPart, ensgs);
                    }).toList();

                    return new TfGroup(tfSymbol, transcriptionFactors, srcDir);
                }).collect(Collectors.toSet());

                Map<String, TfGroup> tfGroupMap = groups.stream().collect(toMap(TfGroup::getSymbol, g -> g));

                logger.trace("Reading log2 fold changes");

                pairingDeseq.forEach((pairing, deseqFile) -> {
                    Map<String, Double> ensgLog2fc;
                    try {
                        ensgLog2fc = readLines(deseqFile).stream().skip(1).map(line -> line.split("\t")).map(
                                split -> Pair.of(split[0].replace("\"", ""), Double.parseDouble(split[1]))).collect(
                                toMap(Pair::getKey, Pair::getValue));
                    } catch (IOException e) {
                        throw new UncheckedIOException(e);
                    }
                    groups.forEach(group -> group.getTranscriptionFactors().forEach(tf -> tf.setLog2fc(pairing,
                            tf.getEnsgs().stream().filter(ensgLog2fc::containsKey).mapToDouble(
                                    ensgLog2fc::get).average().orElse(0.0))));
                });

                logger.trace("Reading mean expression");

                groupMeanExpression.forEach((group, meanFile) -> {
                    Map<String, Double> ensgMeanExpression;
                    try {
                        ensgMeanExpression = readLines(meanFile).stream().skip(1).map(line -> line.split("\t")).map(
                                split -> Pair.of(split[0], Double.parseDouble(split[1]))).collect(
                                toMap(Pair::getKey, Pair::getValue));
                    } catch (IOException e) {
                        throw new UncheckedIOException(e);
                    }
                    groups.forEach(tfGroup -> tfGroup.getTranscriptionFactors().forEach(
                            tf -> tf.setMeanExpression(group,
                                    tf.getEnsgs().stream().filter(ensgMeanExpression::containsKey).mapToDouble(
                                            ensgMeanExpression::get).average().orElse(-1))));
                });

                logger.trace("Reading tpm");

                groupTpm.forEach((group, tpmFile) -> {
                    Map<String, Double> ensgTpm;
                    try {
                        ensgTpm = readLines(tpmFile).stream().skip(1).map(line -> line.split("\t")).map(
                                split -> Pair.of(split[0], Double.parseDouble(split[1]))).collect(
                                toMap(Pair::getKey, Pair::getValue));
                    } catch (IOException e) {
                        throw new UncheckedIOException(e);
                    }
                    groups.forEach(tfGroup -> tfGroup.getTranscriptionFactors().forEach(tf -> tf.setTpm(group,
                            tf.getEnsgs().stream().filter(ensgTpm::containsKey).mapToDouble(
                                    ensgTpm::get).average().orElse(0.0))));
                });

                logger.trace("Reading dcg files");

                hmDcgFiles.forEach((hm, dcgFile) -> {
                    Map<String, Double> symbolDcg;
                    try {
                        symbolDcg = readLines(dcgFile).stream().skip(1).map(line -> line.split("\t")).map(
                                split -> Pair.of(split[0], Double.parseDouble(split[1]))).collect(
                                toMap(Pair::getKey, Pair::getValue));
                    } catch (IOException e) {
                        throw new UncheckedIOException(e);
                    }
                    groups.forEach(group -> group.setRank(hm, symbolDcg.getOrDefault(group.getSymbol(), -1.)));
                });

                File assets = new File(srcDir, "assets");


                Map<String, File> tfBiophysical = Arrays.stream(
                        Objects.requireNonNull(biophysicalLogo.listFiles(file -> file.getName().endsWith(".png")))).map(
                        file -> {
                            String symbol = file.getName().replace(".png", "");
                            return Pair.of(symbol, file);
                        }).collect(toMap(Pair::getKey, Pair::getValue));
                groups.stream().filter(group -> tfBiophysical.containsKey(group.getSymbol())).forEach(group -> {
                    try {
                        group.setBiophysicalLogo(tfBiophysical.get(group.getSymbol()));
                    } catch (IOException e) {
                        throw new UncheckedIOException(e);
                    }
                });

                Map<String, Collection<File>> tfJaspar =
                        Arrays.stream(Objects.requireNonNull(tfSequence.listFiles(File::isDirectory))).map(
                                directory -> {
                                    String symbol = directory.getName();
                                    return Pair.of(symbol,
                                            Arrays.stream(Objects.requireNonNull(directory.listFiles())).filter(
                                                    file -> file.getName().endsWith(".svg")).collect(
                                                    Collectors.toList()));
                                }).collect(toMap(Pair::getKey, Pair::getValue));

                groups.stream().filter(group -> tfJaspar.containsKey(group.getSymbol())).forEach(
                        group -> group.setTfSequence(tfJaspar.get(group.getSymbol())));

                pairingHmHeatmapDir.forEach((pairing, hmHeatmapDir) -> hmHeatmapDir.forEach(
                        (hm, heatmapDir) -> Arrays.stream(Objects.requireNonNull(
                                heatmapDir.listFiles(file -> file.getName().endsWith(".png")))).forEach(file -> {
                            String symbol = file.getName().replace(".png", "");
                            if (tfGroupMap.containsKey(symbol)) {
                                tfGroupMap.get(symbol).setHeatmap(pairing, hm, file);
                            }
                        })));

                pairingHmIgvDir.forEach((pairing, hmIgvDir) -> hmIgvDir.forEach((hm, igvDir) -> Arrays.stream(
                        Objects.requireNonNull(igvDir.listFiles(File::isDirectory))).forEach(directory -> {
                    String symbol = directory.getName();
                    if (tfGroupMap.containsKey(symbol)) {
                        tfGroupMap.get(symbol).setIgv(pairing, hm, directory);
                    }
                })));

                hmDistributionPlots.forEach((hm, distributionPlotDirectory) -> Arrays.stream(Objects.requireNonNull(
                        distributionPlotDirectory.listFiles(f -> f.getName().endsWith(".png")))).forEach(plotFile -> {
                    String symbol = plotFile.getName().replace(".png", "");
                    if (tfGroupMap.containsKey(symbol)) {
                        tfGroupMap.get(symbol).setDistributionPlot(hm, plotFile);
                    }
                }));

                Map<String, Map<String, Map<String, Double>>> hmPairingTfRegressionCoefficient =
                        pairingHmRegressionCoefficients.entrySet().stream().flatMap(pairingHmRegressionCoefficient -> {
                            String pairing = pairingHmRegressionCoefficient.getKey();
                            Map<String, InputFile> hmRegressionCoefficients = pairingHmRegressionCoefficient.getValue();

                            return hmRegressionCoefficients.entrySet().stream().map(hmRegressionCoefficient -> {
                                String hm = hmRegressionCoefficient.getKey();
                                InputFile regressionCoefficient = hmRegressionCoefficient.getValue();
                                Map<String, Double> tfRegressionCoefficient;
                                try {
                                    tfRegressionCoefficient = readLines(regressionCoefficient).stream().skip(1).map(
                                            line -> line.split("\t")).collect(
                                            toMap(split -> split[0], split -> Double.parseDouble(split[1])));
                                } catch (IOException e) {
                                    throw new RuntimeException(e);
                                }

                                return Pair.of(hm, Pair.of(pairing, tfRegressionCoefficient));
                            });
                        }).collect(
                                groupingBy(Pair::getKey, mapping(Pair::getValue, toMap(Pair::getKey, Pair::getValue))));

                final JSONArray coOccurrenceJson;
                if (coOccurrence != null) {
                    try (var reader = new BufferedReader(new FileReader(coOccurrence))) {
                        String[] header = reader.readLine().split("\t");
                        List<JSONObject> objects = reader.lines().map(line -> {
                            String[] split = line.split("\t");
                            String tf = split[0];
                            return Pair.of(tf, Arrays.stream(split).skip(1).mapToDouble(Double::parseDouble).toArray());
                        }).map(pair -> {
                            double[] values = pair.getValue();
                            JSONArray series =
                                    new JSONArray(IntStream.range(0, values.length).mapToObj(i -> new JSONObject() {{
                                        put("name", header[i + 1]);
                                        put("value", values[i]);
                                    }}).collect(toList()));

                            return new JSONObject() {{
                                put("name", pair.getKey());
                                put("series", series);
                            }};
                        }).collect(toList());

                        coOccurrenceJson = new JSONArray(objects);
                    }
                } else {
                    coOccurrenceJson = new JSONArray();
                }

                logger.trace("Reading important loci");
                File importantLociDir = new File(assets, "importantLoci");
                JSONObject importantLociObject = new JSONObject() {{
                    importantLoci.forEach((group, groupInDir) -> {
                        File groupOutDir = new File(importantLociDir, group);
                        JSONObject patternFiles = new JSONObject() {{
                            Arrays.stream(Objects.requireNonNull(groupInDir.listFiles(File::isDirectory))).forEach(
                                    patternInDir -> {
                                        File patternOutDir = new File(groupOutDir, patternInDir.getName());

                                        JSONArray tfFiles = new JSONArray(
                                                Arrays.stream(Objects.requireNonNull(patternInDir.listFiles())).map(
                                                        inFile -> {
                                                            File outFile = new File(patternOutDir, inFile.getName());

                                                            try {
                                                                softLink(outFile, inFile);
                                                            } catch (IOException e) {
                                                                throw new RuntimeException(e);
                                                            }

                                                            return srcDir.toPath().relativize(
                                                                    outFile.toPath()).toString();
                                                        }).toList());

                                        put(patternInDir.getName(), tfFiles);
                                    });
                        }};

                        put(group, patternFiles);
                    });
                }};

                logger.trace("Reading top log2fc");
                File topLog2fcDir = new File(assets, "topLog2fc");
                JSONObject topLog2fcObject = new JSONObject() {{
                    topLog2fc.forEach((pairing, directories) -> {
                        File pairingOutDir = new File(topLog2fcDir, pairing);
                        File upregulatedInDir = directories.getLeft();
                        File downregulatedInDir = directories.getRight();
                        File upregulatedOutDir = new File(pairingOutDir, "upregulated");
                        File downregulatedOutDir = new File(pairingOutDir, "downregulated");
                        put(pairing, new JSONObject() {{
                            put("upregulated", new JSONArray(Arrays.stream(upregulatedInDir.listFiles()).map(inFile -> {
                                File outFile = new File(upregulatedOutDir, inFile.getName());

                                try {
                                    softLink(outFile, inFile);
                                } catch (IOException e) {
                                    throw new RuntimeException(e);
                                }

                                return srcDir.toPath().relativize(outFile.toPath()).toString();
                            }).toList()));

                            put("downregulated",
                                    new JSONArray(Arrays.stream(downregulatedInDir.listFiles()).map(inFile -> {
                                        File outFile = new File(downregulatedOutDir, inFile.getName());

                                        try {
                                            softLink(outFile, inFile);
                                        } catch (IOException e) {
                                            throw new RuntimeException(e);
                                        }

                                        return srcDir.toPath().relativize(outFile.toPath()).toString();
                                    }).toList()));
                        }});
                    });
                }};

                logger.trace("Writing output");

                JSONObject json = new JSONObject() {{
                    put("groups", groups.stream().map(TfGroup::toJSON).toList());
                    put("regressionCoefficients", hmPairingTfRegressionCoefficient);
                    put("configs", configs.getConfigsJSONObject(true, true));
                    put("coOccurrence", coOccurrenceJson);
                    put("importantLoci", importantLociObject);
                    put("topLog2fc", topLog2fcObject);
                }};

                File data = new File(assets, "data.json");

                makeSureFileExists(data);

                try (FileWriter writer = new FileWriter(data)) {
                    json.write(writer, 4, 0);
                }

                return true;
            });
        }};
    }
}
