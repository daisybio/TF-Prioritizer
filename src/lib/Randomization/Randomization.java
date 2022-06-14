package lib.Randomization;

import lib.ExecutableStep;
import tfprio.tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.FileFilters.Filters;
import util.Regions.Region;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;

import static util.FileManagement.extend;
import static util.FileManagement.readLines;

public class Randomization extends ExecutableStep
{
    private final AbstractConfig<File> d_chipAtlas_peakFiles = TFPRIO.configs.chipAtlas.fileStructure.d_peakFiles;
    private final AbstractConfig<File> d_tepic_predictedPeaks =
            TFPRIO.configs.tepic.fileStructure.d_postprocessing_trapPredictedBeds;
    private final AbstractConfig<File> d_input_experimentalPeaks = TFPRIO.configs.igv.pathToTfChipSeq;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>()
        {{
            add(d_chipAtlas_peakFiles);
            add(d_tepic_predictedPeaks);
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
        for (File d_tfGroup : Objects.requireNonNull(d_chipAtlas_peakFiles.get().listFiles(Filters.directoryFilter)))
        {
            String tfGroup = d_tfGroup.getName().replaceAll("\\d+_", "");
            logger.info("Loading data for tf: " + tfGroup);

            Set<Region> chipAtlasRegions = new HashSet<>();
            Set<Region> predictedRegions = new HashSet<>();

            for (File f_bed : Objects.requireNonNull(d_tfGroup.listFiles(Filters.getSuffixFilter(".bed"))))
            {
                chipAtlasRegions.addAll(getRegionsFromBed(f_bed));
            }

            for (File d_group : Objects.requireNonNull(d_tepic_predictedPeaks.get().listFiles(Filters.directoryFilter)))
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
}
