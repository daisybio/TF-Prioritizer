package org.exbio.tfprio.steps.chipAtlas;

import org.apache.commons.lang3.tuple.Pair;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;
import org.exbio.tfprio.lib.Region;
import org.exbio.tfprio.lib.RegionTreeSet;
import org.json.JSONObject;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.exbio.pipejar.util.FileManagement.readLines;

public class ConfusionMatrixes extends ExecutableStep<Configs> {
    public final OutputFile outputFile = addOutput("matrixes");
    private final InputFile chipAtlas;
    private final InputFile chromosomeLengthsFile;
    private final Set<InputFile> predictedRegionDirectories = new HashSet<>();
    private final RequiredConfig<Integer> windowSize = new RequiredConfig<>(configs.tepic.windowSize);
    private final int margin;

    public ConfusionMatrixes(Configs configs, OutputFile chipAtlas, OutputFile chromosomeLengthsFile,
                             Map<String, Map<String, OutputFile>> predictedRegionDirectories) {
        super(configs, false, chipAtlas, chromosomeLengthsFile);

        this.chipAtlas = addInput(chipAtlas);
        this.chromosomeLengthsFile = addInput(chromosomeLengthsFile);

        OutputFile predictedDir = new OutputFile(inputDirectory, "predicted");
        predictedRegionDirectories.forEach((group, hmMap) -> {
            OutputFile groupDir = new OutputFile(predictedDir, group);
            hmMap.forEach((hm, hmInDir) -> this.predictedRegionDirectories.add(addInput(groupDir, hmInDir)));
        });

        margin = windowSize.get() / 2;
    }

    private static Collection<Region> getRandomEntries(Collection<Region> regions, int numberOfEntries) {
        List<Region> allRegions = new ArrayList<>(regions);

        int numberOfAvailableRegions = regions.size();

        if (numberOfEntries < numberOfAvailableRegions) {
            // Prevents duplicates
            Collections.shuffle(allRegions);
            return new ArrayList<>(allRegions.subList(0, numberOfEntries));
        } else {
            // Duplicates are allowed, otherwise there would be too few entries
            return IntStream.range(0, numberOfEntries).mapToObj(
                    i -> allRegions.get(new Random().nextInt(allRegions.size()))).toList();
        }
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            logger.trace("Reading chromosome lengths");
            Map<String, Integer> chromosomeLengths;
            try {
                chromosomeLengths = readLines(chromosomeLengthsFile).stream().map(line -> line.split("\t")).collect(
                        Collectors.toMap(line -> line[0], line -> Integer.parseInt(line[1])));
            } catch (IOException e) {
                throw new RuntimeException(e);
            }


            Arrays.stream(Objects.requireNonNull(chipAtlas.listFiles(file -> file.getName().endsWith(".bed")))).forEach(
                    bedFile -> add(() -> {
                        String tf = bedFile.getName().split("_")[0];
                        final RegionTreeSet chipRegions;
                        final RegionTreeSet predictedRegions;

                        try {
                            chipRegions =
                                    new RegionTreeSet(readLines(bedFile).stream().skip(1).map(Region::new).toList());
                        } catch (IOException e) {
                            throw new UncheckedIOException(e);
                        }

                        predictedRegions = new RegionTreeSet(predictedRegionDirectories.stream().map(
                                directory -> new File(directory, tf + ".bed")).filter(File::exists).flatMap(file -> {
                            try {
                                return readLines(file).stream().skip(1).map(Region::new);
                            } catch (IOException e) {
                                throw new RuntimeException(e);
                            }
                        }).toList());

                        int truePositives = (int) predictedRegions.stream().filter(
                                predictedRegion -> chipRegions.hasOverlap(predictedRegion, margin)).count();
                        int falsePositives = predictedRegions.size() - truePositives;

                        final RegionTreeSet notPredictedRegions = new RegionTreeSet(predictedRegions.stream().map(
                                predictedRegion -> Pair.of(predictedRegions.lower(predictedRegion),
                                        predictedRegion)).filter(pair -> pair.getLeft() != null).filter(
                                pair -> pair.getLeft().getChromosome().equals(pair.getRight().getChromosome())).map(
                                pair -> new Region(pair.getLeft().getChromosome(), pair.getLeft().getEnd() + 1,
                                        pair.getRight().getStart() - 1)).toList());

                        chromosomeLengths.entrySet().stream().map(entry -> Pair.of(new Region(entry.getKey(), 0, 0),
                                new Region(entry.getKey(), entry.getValue(), entry.getValue()))).flatMap(pair -> {
                            Region first = predictedRegions.ceiling(pair.getLeft());
                            Region last = predictedRegions.floor(pair.getRight());

                            Set<Region> regions = new HashSet<>();

                            if (first != null) {
                                regions.add(new Region(first.getChromosome(), 0, first.getStart() - 1));
                            }
                            if (last != null) {
                                regions.add(new Region(last.getChromosome(), last.getEnd() + 1,
                                        chromosomeLengths.get(last.getChromosome())));
                            }

                            return regions.stream();
                        }).forEach(notPredictedRegions::add);

                        final RegionTreeSet notPredictedRegionsSplitted =
                                new RegionTreeSet(notPredictedRegions.stream().flatMap(region -> {
                                    int start = region.getStart();
                                    int end = region.getEnd();
                                    int length = end - start + 1;
                                    int numberOfWindows = (length + windowSize.get() - 1) / windowSize.get();

                                    return IntStream.range(0, numberOfWindows).mapToObj(
                                            i -> new Region(region.getChromosome(), start + i * windowSize.get(),
                                                    Math.min(start + (i + 1) * windowSize.get() - 1, end)));
                                }).toList());

                        final Collection<Region> notPredictedRegionsRandom =
                                getRandomEntries(notPredictedRegionsSplitted, predictedRegions.size());

                        int falseNegatives = (int) notPredictedRegionsRandom.stream().filter(
                                predictedRegion -> chipRegions.hasOverlap(predictedRegion, margin)).count();
                        int trueNegatives = notPredictedRegionsRandom.size() - falseNegatives;

                        File matrixFile = new File(outputFile, tf + ".json");
                        JSONObject matrix = new JSONObject() {{
                            put("truePositives", truePositives);
                            put("falsePositives", falsePositives);
                            put("falseNegatives", falseNegatives);
                            put("trueNegatives", trueNegatives);
                        }};

                        try (var writer = new FileWriter(matrixFile)) {
                            matrix.write(writer, 4, 0);
                        } catch (IOException e) {
                            throw new UncheckedIOException(e);
                        }

                        return true;
                    }));
        }};
    }
}
