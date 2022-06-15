package lib.Randomization;

import lib.ExecutableStep;
import org.json.JSONObject;
import tfprio.tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.FileFilters.Filters;
import util.RegionSearchTree.ChromosomeRegionTrees;
import util.Regions.Region;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import static util.FileManagement.*;

public class Randomization extends ExecutableStep
{
    private final AbstractConfig<File> f_input_chromosomeLengths =
            TFPRIO.configs.deSeq2.fileStructure.f_preprocessing_chromosomeLengths;
    private final AbstractConfig<File> d_chipAtlas_peakFiles = TFPRIO.configs.chipAtlas.fileStructure.d_peakFiles;
    private final AbstractConfig<File> d_tepic_predictedPeaks =
            TFPRIO.configs.tepic.fileStructure.d_postprocessing_trapPredictedBeds;
    private final AbstractConfig<File> d_input_experimentalPeaks = TFPRIO.configs.igv.pathToTfChipSeq;
    private final AbstractConfig<Integer> windowSize = TFPRIO.configs.tepic.windowSize;

    private final Map<String, Integer> chromosomeLengths = new HashMap<>();

    private static List<Region> getRandomEntries(List<Region> regions, int numberOfEntries)
    {
        List<Region> selectedRegions;
        int numberOfAvailableRegions = regions.size();
        if (numberOfEntries < numberOfAvailableRegions)
        {
            Collections.shuffle(regions);
            selectedRegions = new ArrayList<>(regions.subList(0, numberOfEntries));
        } else
        {
            selectedRegions = new ArrayList<>();
            Random random = new Random();

            for (int i = 0; i < numberOfEntries; i++)
            {
                selectedRegions.add(regions.get(random.nextInt(numberOfAvailableRegions)));
            }
        }

        return selectedRegions;
    }

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>()
        {{
            add(d_chipAtlas_peakFiles);
            add(d_tepic_predictedPeaks);
            add(f_input_chromosomeLengths);
        }};
    }

    @Override public Set<GeneratedFileStructure> getCreatedFileStructure()
    {
        return new HashSet<>();
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>();
    }

    @Override protected Set<AbstractConfig<?>> getOptionalConfigs()
    {
        return new HashSet<>()
        {{
            add(d_input_experimentalPeaks);
        }};
    }

    @Override protected void execute()
    {
        try
        {
            JSONObject json = new JSONObject(readFile(f_input_chromosomeLengths.get()));
            json.keySet().forEach(key -> chromosomeLengths.put(key, json.getInt(key)));
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        System.out.println(chromosomeLengths);

        for (File d_tfGroup : Objects.requireNonNull(d_chipAtlas_peakFiles.get().listFiles(Filters.directoryFilter)))
        {
            {
                String tfGroup = d_tfGroup.getName().replaceAll("\\d+_", "");
                logger.info("Loading data for tf: " + tfGroup);

                Set<Region> chipAtlasRegions = new HashSet<>();
                Set<Region> predictedRegions = new HashSet<>();

                for (File f_bed : Objects.requireNonNull(d_tfGroup.listFiles(Filters.getSuffixFilter(".bed"))))
                {
                    chipAtlasRegions.addAll(getRegionsFromBed(f_bed));
                }

                for (File d_group : Objects.requireNonNull(
                        d_tepic_predictedPeaks.get().listFiles(Filters.directoryFilter)))
                {
                    for (File d_hm : Objects.requireNonNull(d_group.listFiles(Filters.directoryFilter)))
                    {
                        File f_bed;
                        if (tfGroup.contains(".."))
                        {
                            f_bed = extend(d_hm, tfGroup.replace("..", "::") + "_merged.bed");
                        } else
                        {
                            f_bed = extend(d_hm, tfGroup + ".bed");
                        }
                        if (f_bed.exists())
                        {
                            predictedRegions.addAll(getRegionsFromBed(f_bed));
                        }
                    }
                }

                if (predictedRegions.size() == 0)
                {
                    logger.warn("No predicted data found for tf: " + tfGroup);
                } else
                {
                    ConfusionMatrix chipVsPredicted = getConfusionMatrix(chipAtlasRegions, predictedRegions);
                    logger.info("ChipAtlas vs predicted: " + chipVsPredicted);
                }

                if (d_input_experimentalPeaks.isSet())
                {
                    Set<Region> experimentalRegions = new HashSet<>();

                    for (File d_group : Objects.requireNonNull(
                            d_input_experimentalPeaks.get().listFiles(Filters.directoryFilter)))
                    {
                        File d_tf = extend(d_group, tfGroup);
                        if (!d_tf.exists())
                        {
                            continue;
                        }
                        for (File f_peaks : Objects.requireNonNull(d_tf.listFiles(Filters.fileFilter)))
                        {
                            experimentalRegions.addAll(getRegionsFromBed(f_peaks));
                        }
                    }

                    if (experimentalRegions.size() == 0)
                    {
                        logger.warn("No experimental data found for tf: " + tfGroup);
                    } else
                    {
                        ConfusionMatrix experimentalVsPredicted =
                                getConfusionMatrix(experimentalRegions, predictedRegions);
                        logger.info("Experimental vs predicted: " + experimentalVsPredicted);

                        Set<Region> combined = new HashSet<>()
                        {{
                            addAll(experimentalRegions);
                            addAll(chipAtlasRegions);
                        }};
                        ConfusionMatrix combinedVsPredicted = getConfusionMatrix(combined, predictedRegions);
                        logger.info("Combined vs predicted: " + combinedVsPredicted);
                    }
                }
            }
        }
    }

    private Set<Region> getRegionsFromBed(File f_bed)
    {
        Set<Region> regions = new HashSet<>();
        try
        {
            List<String> lines = readLines(f_bed);
            lines.stream().filter(line -> !line.startsWith("track")).map(line -> line.split("\t")).map(Region::new)
                    .forEach(regions::add);
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        return regions;
    }

    private ConfusionMatrix getConfusionMatrix(Iterable<Region> groundTruth, Iterable<Region> comparison)
    {
        ChromosomeRegionTrees groundTruthTrees = new ChromosomeRegionTrees();
        groundTruthTrees.addAllOptimized(groundTruth, true);
        int margin = windowSize.get() / 2;

        Function<Region, Region> getRegionWithMargin =
                region -> new Region(region.getChromosome(), region.getStart() - margin, region.getEnd() + margin);

        int tp = 0, fp = 0, fn = 0, tn = 0;

        for (Region entry : comparison)
        {
            Region searchRegion = getRegionWithMargin.apply(entry);
            if (groundTruthTrees.hasOverlap(searchRegion))
            {
                tp++;
            } else
            {
                fp++;
            }
        }

        ChromosomeRegionTrees comparisonTrees = new ChromosomeRegionTrees();
        comparisonTrees.addAllOptimized(comparison, true);

        for (Region entry : groundTruth)
        {
            Region searchRegion = getRegionWithMargin.apply(entry);
            if (!comparisonTrees.hasOverlap(searchRegion))
            {
                fn++;
            }
        }

        Set<Region> comparisonRegionsWithWindow = new HashSet<>()
        {{
            addAll(StreamSupport.stream(comparison.spliterator(), false).map(getRegionWithMargin)
                    .collect(Collectors.toSet()));
        }};

        ChromosomeRegionTrees comparisonWithWindows = new ChromosomeRegionTrees();
        comparisonWithWindows.addAllOptimized(comparisonRegionsWithWindow, true);
        Map<String, List<Region>> sortedComparisonsWithWindow = comparisonWithWindows.getAllRegionsSorted();
        List<Region> notComparedWindows = new ArrayList<>();

        sortedComparisonsWithWindow.forEach((chromosome, regions) ->
        {
            List<Region> emptyRegions = new ArrayList<>();
            int nextStart = 0;

            for (Region region : regions)
            {
                emptyRegions.add(new Region(chromosome, nextStart, region.getStart() - 1));
                nextStart = region.getEnd() + 1;
            }

            emptyRegions.forEach(region ->
            {
                int length = region.getEnd() - region.getStart();

                for (int end = windowSize.get(); end < length + windowSize.get(); end += windowSize.get() + 1)
                {
                    if (end > length)
                    {
                        end = length;
                    }

                    int start = Math.max(end - windowSize.get(), 0);
                    notComparedWindows.add(new Region(chromosome, start, end));
                }
            });
        });

        List<Region> randomWindows = getRandomEntries(notComparedWindows,
                (int) StreamSupport.stream(comparison.spliterator(), false).count());

        for (Region randomWindow : randomWindows)
        {
            if (groundTruthTrees.hasOverlap(randomWindow))
            {
                fn++;
            } else
            {
                tn++;
            }
        }

        ConfusionMatrix matrix = new ConfusionMatrix();
        matrix.setTp(tp);
        matrix.setFp(fp);
        matrix.setFn(fn);
        matrix.setTn(tn);

        return matrix;
    }
}
