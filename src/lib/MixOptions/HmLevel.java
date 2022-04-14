package lib.MixOptions;

import lib.ExecutableStep;
import tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;

import java.io.File;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class HmLevel extends ExecutableStep
{
    private final AbstractConfig<File> d_input = TFPRIO.configs.mixOptions.fileStructure.d_preprocessingHmMix;
    private final GeneratedFileStructure d_output = TFPRIO.configs.mixOptions.fileStructure.d_hmMix;
    private final AbstractConfig<String> option = TFPRIO.configs.mixOptions.option;
    private final AbstractConfig<Integer> occurrenceIntersection = TFPRIO.configs.mixOptions.occurrenceIntersection;

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
        StaticMethods.runMixOption(d_input.get(), d_output.get(), "HM_LEVEL", logger, executorService, option,
                occurrenceIntersection);
    }
}
