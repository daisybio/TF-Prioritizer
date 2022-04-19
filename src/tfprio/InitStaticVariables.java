package tfprio;

import lib.ExecutableStep;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.FileFilters.Filters;
import util.MapSymbolAndEnsg;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Objects;
import java.util.Set;

import static util.FileManagement.findValueInTable;

public class InitStaticVariables extends ExecutableStep
{
    // Configs for map
    private final GeneratedFileStructure f_map = TFPRIO.configs.deSeq2.fileStructure.f_mapping;
    private final AbstractConfig<File> f_scriptTemplate = TFPRIO.configs.scriptTemplates.f_mapping;
    private final GeneratedFileStructure f_script = TFPRIO.configs.deSeq2.fileStructure.f_mappingScript;
    private final AbstractConfig<File> f_geneIDs = TFPRIO.configs.deSeq2.inputGeneID;
    private final AbstractConfig<File> d_samples = TFPRIO.configs.deSeq2.inputDirectory;
    private final AbstractConfig<File> f_batches = TFPRIO.configs.deSeq2.batchFile;
    private final AbstractConfig<String> datasetSpecies = TFPRIO.configs.deSeq2.biomartDatasetSpecies;
    private final AbstractConfig<String> datasetSymbolColumn = TFPRIO.configs.deSeq2.biomartDatasetSymbolColumn;

    // Configs for inputDirectory
    private final AbstractConfig<File> inputDirectory = TFPRIO.configs.tepic.inputDirectory;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(Arrays.asList(f_scriptTemplate, f_geneIDs, inputDirectory, d_samples))
        {{
            if (f_batches.isSet())
            {
                add(f_batches);
            }
        }};
    }

    @Override public Set<GeneratedFileStructure> getCreatedFileStructure()
    {
        return new HashSet<>(Arrays.asList(f_map, f_script));
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>(Arrays.asList(datasetSpecies, datasetSymbolColumn));
    }

    @Override protected void execute()
    {
        TFPRIO.mapSymbolAndEnsg = new MapSymbolAndEnsg();
        initExistingGroups();
        initSampleGroupMap();
        initSampleBatchMap();
    }

    private void initSampleGroupMap()
    {
        for (File d_group : Objects.requireNonNull(d_samples.get().listFiles(Filters.directoryFilter)))
        {
            for (File f_sample : Objects.requireNonNull(d_group.listFiles(Filters.fileFilter)))
            {
                String sampleName = f_sample.getName().substring(0, f_sample.getName().lastIndexOf("."));
                TFPRIO.sample_group.put(sampleName, d_group.getName());
            }
        }

        logger.info("Samples to groups: " + TFPRIO.sample_group);
    }

    private void initSampleBatchMap()
    {
        if (f_batches.isSet())
        {
            boolean foundAll = true;

            for (String sample : TFPRIO.sample_group.keySet())
            {
                try
                {
                    String batch;
                    batch = findValueInTable(sample, 0, 1, f_batches.get(), "\t", true);
                    TFPRIO.sample_batch.put(sample, batch);
                } catch (FileNotFoundException e)
                {
                    logger.error(e.getMessage());
                } catch (NoSuchFieldException e)
                {
                    foundAll = false;
                    logger.warn(e.getMessage());
                }
            }

            if (!foundAll)
            {
                logger.error(
                        "Could not find all samples in the batch file. Fore more information inspect the previous " +
                                "warnings.");
            } else
            {
                logger.info("Sample to batch: " + TFPRIO.sample_batch);
            }
        } else
        {
            TFPRIO.sample_batch = null;
        }
    }

    private void initExistingGroups()
    {
        for (File d_group : Objects.requireNonNull(inputDirectory.get().listFiles(Filters.directoryFilter)))
        {
            TFPRIO.groupsToHms.put(d_group.getName(), new HashSet<>());
            for (File d_hm : Objects.requireNonNull(d_group.listFiles(Filters.directoryFilter)))
            {
                TFPRIO.existingHms.add(d_hm.getName());
                TFPRIO.groupsToHms.get(d_group.getName()).add(d_hm.getName());
            }
        }

        for (String first : TFPRIO.groupsToHms.keySet())
        {
            for (String second : TFPRIO.groupsToHms.keySet())
            {
                if (first.compareTo(second) >= 0)
                {
                    continue;
                }

                String combination = first + "_" + second;

                TFPRIO.groupCombinationsToHms.put(combination, new HashSet<>());

                for (String hm : TFPRIO.groupsToHms.get(first))
                {
                    if (TFPRIO.groupsToHms.get(second).contains(hm))
                    {
                        TFPRIO.groupCombinationsToHms.get(combination).add(hm);
                    }
                }
            }
        }

        logger.info("Existing hms: " + TFPRIO.existingHms);
        logger.info("Groups to hms: " + TFPRIO.groupsToHms);
        logger.info("Group combinations to hms: " + TFPRIO.groupCombinationsToHms);
    }

    @Override protected boolean mayBeSkipped()
    {
        return false;
    }
}
