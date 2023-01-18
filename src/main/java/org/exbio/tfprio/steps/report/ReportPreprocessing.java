package org.exbio.tfprio.steps.report;

import org.apache.commons.lang3.tuple.Pair;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.json.JSONObject;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import static java.util.stream.Collectors.toMap;
import static org.exbio.pipejar.util.FileManagement.*;

public class ReportPreprocessing extends ExecutableStep {
    public final OutputFile outputFile = addOutput("report");
    private final InputFile ensgSymbol;
    private final InputFile angularInput;
    private final Map<String, InputFile> hmDcgFiles = new HashMap<>();
    private final Map<String, InputFile> pairingDeseq = new HashMap<>();
    private final Map<String, InputFile> groupMeanExpression = new HashMap<>();
    private final Map<String, InputFile> groupTpm = new HashMap<>();
    private final Map<String, Map<String, InputFile>> pairingHmHeatmapDir = new HashMap<>();
    private final InputFile biophysicalLogo;
    private final InputFile tfSequence;

    public ReportPreprocessing(OutputFile ensgSymbol, Map<String, OutputFile> hmDcgFiles,
                               Map<String, OutputFile> deseqResults, Map<String, OutputFile> groupMeanExpression,
                               Map<String, OutputFile> groupTpm, OutputFile biophysicalLogo, OutputFile tfSequence,
                               Map<String, Map<String, OutputFile>> pairingHmHeatmapDir) {
        super(false, hmDcgFiles.values(), deseqResults.values(), groupMeanExpression.values(), groupTpm.values(),
                List.of(ensgSymbol),
                pairingHmHeatmapDir.values().stream().flatMap(m -> m.values().stream()).collect(Collectors.toList()));

        this.ensgSymbol = addInput(ensgSymbol);
        this.biophysicalLogo = addInput(biophysicalLogo);
        this.tfSequence = addInput(tfSequence);
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
                                split -> split.length > 1).collect(Collectors.groupingBy(split -> split[1],
                                Collectors.mapping(split -> split[0], Collectors.toList())));

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


                logger.trace("Writing output");

                JSONObject json = new JSONObject() {{
                    put("groups", groups.stream().map(TfGroup::toJSON).toList());
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
