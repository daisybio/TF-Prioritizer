package lib.MixOptions;

import lib.ExecutableStep;
import tfprio.TFPRIO;
import util.Configs.Config;

import java.io.File;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class SampleLevel extends ExecutableStep
{
    private final Config<File> d_input = TFPRIO.configs.mixOptions.fileStructure.d_sampleMixPreprocessing;
    private final Config<File> d_output = TFPRIO.configs.mixOptions.fileStructure.d_sampleMix;
    private final Config<String> option = TFPRIO.configs.mixOptions.option;
    private final Config<Integer> occurrenceIntersection = TFPRIO.configs.mixOptions.occurrenceIntersection;

    @Override protected Set<Config<File>> getRequiredFileStructure()
    {
        return new HashSet<>(List.of(d_input));
    }

    @Override protected Set<Config<File>> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(d_output));
    }

    @Override protected Set<Config<?>> getRequiredConfigs()
    {
        return new HashSet<>(Arrays.asList(option, occurrenceIntersection));
    }

    @Override protected void updateInputDirectory()
    {
        TFPRIO.latestInputDirectory = d_output;
    }

    @Override protected void execute()
    {
        StaticMethods.runMixOption(d_input.get(), d_output.get(), "SAMPLE_LEVEL", logger, executorService, option,
                occurrenceIntersection);
    }
}