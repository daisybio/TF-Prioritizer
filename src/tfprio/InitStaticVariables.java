package tfprio;

import lib.ExecutableStep;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.FileFilters.Filters;
import util.MapSymbolAndEnsg;

import java.io.File;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Objects;
import java.util.Set;

public class InitStaticVariables extends ExecutableStep
{
    // Configs for map
    private final GeneratedFileStructure f_map = TFPRIO.configs.deSeq2.fileStructure.f_mapping;
    private final AbstractConfig<File> f_scriptTemplate = TFPRIO.configs.scriptTemplates.f_mapping;
    private final GeneratedFileStructure f_script = TFPRIO.configs.deSeq2.fileStructure.f_mappingScript;
    private final AbstractConfig<File> f_geneIDs = TFPRIO.configs.deSeq2.inputGeneID;
    private final AbstractConfig<String> datasetSpecies = TFPRIO.configs.deSeq2.biomartDatasetSpecies;
    private final AbstractConfig<String> datasetSymbolColumn = TFPRIO.configs.deSeq2.biomartDatasetSymbolColumn;

    // Configs for inputDirectory
    private final AbstractConfig<File> inputDirectory = TFPRIO.configs.tepic.inputDirectory;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(Arrays.asList(f_scriptTemplate, f_geneIDs, inputDirectory));
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
