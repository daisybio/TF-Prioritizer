package lib.MixOptions;

import lib.ExecutableStep;
import tfprio.tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.FileFilters.Filters;

import java.io.File;
import java.util.*;

import static util.FileManagement.extend;

public class HmLevelPreprocess extends ExecutableStep
{
    private final AbstractConfig<File> d_input = TFPRIO.configs.mixOptions.fileStructure.d_sampleMix;
    private final GeneratedFileStructure d_output = TFPRIO.configs.mixOptions.fileStructure.d_preprocessingHmMix;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(List.of(d_input));
    }

    @Override public Set<GeneratedFileStructure> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(d_output));
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>();
    }

    @Override protected void execute()
    {
        logger.info("Preprocess sample unions for HM mix");
        logger.info("Identify possible groups with same histone modifications");

        HashMap<String, ArrayList<String>> groups_histoneModifications = new HashMap<>();
        HashMap<String, ArrayList<String>> deleted_groups = new HashMap<>();
        HashSet<String> available_histoneModifications = new HashSet<>();

        //identify possible timepoints
        for (File d_group : Objects.requireNonNull(d_input.get().listFiles(Filters.directoryFilter)))
        {
            String group = d_group.getName();
            ArrayList<String> hmList = new ArrayList<>();

            for (File d_hm : Objects.requireNonNull(d_group.listFiles(Filters.directoryFilter)))
            {
                hmList.add(d_hm.getName());
                available_histoneModifications.add(d_hm.getName());
            }
            groups_histoneModifications.put(group, hmList);
        }


        for (String group : groups_histoneModifications.keySet())
        {
            if (groups_histoneModifications.get(group).size() < available_histoneModifications.size())
            {
                deleted_groups.put(group, groups_histoneModifications.get(group));
                continue;
            }

            if (!groups_histoneModifications.get(group).containsAll(available_histoneModifications))
            {
                deleted_groups.put(group, groups_histoneModifications.get(group));
            }
        }

        deleted_groups.keySet().forEach(groups_histoneModifications::remove);

        StringBuilder sb_foundMixingGroups = new StringBuilder();

        sb_foundMixingGroups.append("Can perform complete mix for HMs (");

        available_histoneModifications.forEach(
                histoneModification -> sb_foundMixingGroups.append(histoneModification).append(" "));

        sb_foundMixingGroups.append(") in timepoints (");

        groups_histoneModifications.keySet().forEach(group -> sb_foundMixingGroups.append(group).append(" "));

        sb_foundMixingGroups.append("). Can perform part-mix or no-mix for timepoints (");

        deleted_groups.keySet().forEach(deletedGroup -> sb_foundMixingGroups.append(deletedGroup).append(" "));

        sb_foundMixingGroups.append(").");
        logger.info(sb_foundMixingGroups.toString());


        for (File d_group : Objects.requireNonNull(d_input.get().listFiles(Filters.directoryFilter)))
        {
            String group = d_group.getName();
            File d_output_hmPreprocessing = extend(d_output.get(), group);

            for (File d_hm : Objects.requireNonNull(d_group.listFiles(Filters.directoryFilter)))
            {
                File d_output_hmPreprocessingMix = extend(d_output_hmPreprocessing, "MIX");

                for (File f_sample : Objects.requireNonNull(d_hm.listFiles(Filters.fileFilter)))
                {
                    StaticMethods.splitFileByChromosome(f_sample, d_output_hmPreprocessingMix, logger);
                }
            }
        }
    }
}
