package lib.MixOptions;

import lib.ExecutableStep;
import tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.FileFilters.Filters;

import java.io.*;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;

import static util.FileManagement.extend;

public class SampleLevelPreprocess extends ExecutableStep
{
    private final AbstractConfig<File> d_input = TFPRIO.configs.mixOptions.fileStructure.d_preprocessingCheckChr;
    private final AbstractConfig<File> d_output = TFPRIO.configs.mixOptions.fileStructure.d_sampleMixPreprocessing;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(List.of(d_input));
    }

    @Override public Set<AbstractConfig<File>> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(d_output));
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>();
    }

    @Override protected void execute()
    {
        for (File d_group : Objects.requireNonNull(d_input.get().listFiles(Filters.directoryFilter)))
        {
            String group = d_group.getName();
            File d_outputGroup = extend(d_output.get(), group);

            for (File d_hm : Objects.requireNonNull(d_group.listFiles(Filters.directoryFilter)))
            {
                File d_outputGroupHm = extend(d_outputGroup, d_hm.getName());

                for (File f_sample : Objects.requireNonNull(d_hm.listFiles(Filters.fileFilter)))
                {
                    executorService.execute(() ->
                    {
                        StaticMethods.splitFileByChromosome(f_sample, d_outputGroupHm, logger);
                    });
                }
            }
        }
    }
}
