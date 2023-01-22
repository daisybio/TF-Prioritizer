package org.exbio.tfprio.steps.chipAtlas;

import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.net.URL;
import java.nio.file.Files;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.exbio.pipejar.util.FileManagement.*;

public class GetData extends ExecutableStep<Configs> {
    public final OutputFile outputFile = addOutput("chipAtlas");
    private final InputFile dataList;
    private final InputFile dcgFile;
    private final RequiredConfig<List<String>> tissueTypes = new RequiredConfig<>(configs.chipAtlas.tissueTypes);
    private final RequiredConfig<String> genomeVersionColName =
            new RequiredConfig<>(configs.chipAtlas.genomeVersionColName);
    private final RequiredConfig<String> antigenClassColName =
            new RequiredConfig<>(configs.chipAtlas.antigenClassColName);
    private final RequiredConfig<String> antigenColName = new RequiredConfig<>(configs.chipAtlas.antigenRegexColName);
    private final RequiredConfig<String> cellTypeClassColName =
            new RequiredConfig<>(configs.chipAtlas.cellTypeClassColName);
    private final RequiredConfig<String> urlColName = new RequiredConfig<>(configs.chipAtlas.urlColName);
    private final RequiredConfig<String> thresholdColName = new RequiredConfig<>(configs.chipAtlas.thresholdColName);
    private final RequiredConfig<String> genomeVersion = new RequiredConfig<>(configs.chipAtlas.genomeVersion);

    public GetData(Configs configs, OutputFile dataList, OutputFile dcgFile) {
        super(configs, false, dataList);

        this.dataList = addInput(dataList);
        this.dcgFile = addInput(dcgFile);
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            Set<RequiredConfig<String>> colNames =
                    Set.of(genomeVersionColName, antigenClassColName, antigenColName, cellTypeClassColName, urlColName,
                            thresholdColName);

            Set<String> tfNames;
            Map<String, String> tfUrl;

            try (var reader = new BufferedReader(new FileReader(dataList))) {
                String[] header = reader.readLine().split(",");

                Map<RequiredConfig<String>, Integer> colNameIndex = colNames.stream().map(colNameConfigs -> {
                    int index = IntStream.range(0, header.length).filter(
                            i -> header[i].equals(colNameConfigs.get())).findFirst().orElseThrow();
                    return new AbstractMap.SimpleEntry<>(colNameConfigs, index);
                }).collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));


                if (colNameIndex.size() != colNames.size()) {
                    throw new RuntimeException("Not all column names found in header. Missing: " +
                            colNames.stream().filter(regex -> !colNameIndex.containsKey(regex)).map(
                                    RequiredConfig::get).collect(Collectors.joining(", ")));
                }

                tfNames = readLines(dcgFile).stream().skip(1).map(line -> line.split("\t")[0]).collect(
                        Collectors.toSet());

                logger.trace("Starting to read data list");

                // Using the sophisticated regex right from the start leads to performance issues
                tfUrl = reader.lines().map(line -> Pair.of(line, line.split(","))).filter(
                        pair -> pair.getValue()[colNameIndex.get(genomeVersionColName)].equalsIgnoreCase(
                                genomeVersion.get())).filter(
                        pair -> pair.getValue()[colNameIndex.get(antigenClassColName)].equals("TFs and others")).filter(
                        split -> tissueTypes.get().stream().anyMatch(tissueType -> tissueType.equalsIgnoreCase(
                                split.getValue()[colNameIndex.get(cellTypeClassColName)]))).map(pair -> {
                    String tf = pair.getValue()[colNameIndex.get(antigenColName)].toUpperCase();
                    int threshold = Integer.parseInt(pair.getValue()[colNameIndex.get(thresholdColName)]);
                    String[] trueSplit = pair.getKey().split(",(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)");
                    String url = trueSplit[colNameIndex.get(urlColName)];
                    return Pair.of(tf, Pair.of(threshold, url));
                }).collect(Collectors.groupingBy(Pair::getKey,
                        Collectors.minBy(Comparator.comparing(p -> p.getValue().getKey())))).entrySet().stream().filter(
                        entry -> entry.getValue().isPresent()).collect(
                        Collectors.toMap(Map.Entry::getKey, e -> e.getValue().get().getValue().getValue()));

                if (tfUrl.isEmpty()) {
                    logger.warn("No TFs found for given tissue types");
                } else {
                    logger.debug("Found {} TFs for given tissue types", tfUrl.size());
                }
            } catch (IOException e) {
                throw new RuntimeException(e);
            }

            Map<String, Set<String>> tfSplits =
                    tfNames.stream().map(tf -> Pair.of(tf, new HashSet<>(Arrays.asList(tf.split("::"))))).peek(
                            pair -> pair.getValue().retainAll(tfUrl.keySet())).filter(
                            pair -> !pair.getValue().isEmpty()).collect(Collectors.toMap(Pair::getKey, Pair::getValue));

            logger.debug("Found {} TFs which are available from both ChipAtlas and DCG file", tfSplits.size());

            Set<Pair<String, String>> splitUrls = tfSplits.values().stream().flatMap(Collection::stream).distinct().map(
                    tf -> Pair.of(tf, tfUrl.get(tf))).collect(Collectors.toSet());

            splitUrls.forEach(splitUrl -> add(() -> {
                String tf = splitUrl.getKey();
                String url = splitUrl.getValue();
                File bedFile = new File(outputFile, tf + "_" + url.substring(url.lastIndexOf('/') + 1));
                int attempt = 1;

                while (!bedFile.exists()) {
                    logger.debug("Attempt " + attempt + " to download chip atlas bed file");
                    try {
                        makeSureFileExists(bedFile);
                        IOUtils.copy(new URL(url), bedFile);

                        long size = Files.size(bedFile.toPath());

                        if (size < 300) {
                            logger.warn("File {} is too small ({} bytes)", bedFile, size);
                            throw new IOException("File too small");
                        }
                    } catch (IOException e) {
                        deleteFileStructure(bedFile);
                        attempt++;
                        if (attempt > 3) {
                            throw new RuntimeException("Could not download file " + url);
                        }
                    }
                }

                return true;
            }));
        }};
    }
}
