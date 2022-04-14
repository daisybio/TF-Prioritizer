package lib.Tgene;

import lib.ExecutableStep;
import tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.FileFilters.Filters;

import java.io.*;
import java.util.*;

import static util.FileManagement.*;
import static util.ScriptExecution.executeAndWait;

public class RunTgene extends ExecutableStep
{
    private AbstractConfig<File> d_input;
    private final AbstractConfig<File> f_input_gtf = TFPRIO.configs.tgene.fileStructure.f_transcripts_gtf;

    private final AbstractConfig<File> d_output = TFPRIO.configs.tgene.fileStructure.d_output;

    private final AbstractConfig<File> pathToTgeneExecutable = TFPRIO.configs.tgene.pathToExecutable;

    private final AbstractConfig<Boolean> noClosestLocus = TFPRIO.configs.tgene.noClosestLocus;
    private final AbstractConfig<Boolean> noClosestTss = TFPRIO.configs.tgene.noClosestTss;
    private final AbstractConfig<Integer> maxLinkDistance = TFPRIO.configs.tgene.maxLinkDistance;
    private final AbstractConfig<Double> pValue = TFPRIO.configs.tgene.pValue;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(Arrays.asList(d_input, pathToTgeneExecutable, f_input_gtf));
    }

    @Override protected Set<AbstractConfig<File>> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(d_output));
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>(Arrays.asList(noClosestLocus, noClosestTss, maxLinkDistance, pValue));
    }

    @Override protected void updateInputDirectory()
    {
        d_input = TFPRIO.latestInputDirectory;
    }

    @Override protected void execute()
    {
        File tgeneExecutable = extend(pathToTgeneExecutable.get(), "bin", "tgene");

        for (File d_group : Objects.requireNonNull(d_input.get().listFiles(Filters.directoryFilter)))
        {
            for (File d_hm : Objects.requireNonNull(d_group.listFiles(Filters.directoryFilter)))
            {
                for (File f_sample : Objects.requireNonNull(d_hm.listFiles(Filters.fileFilter)))
                {
                    try
                    {
                        // Prevent crashes due to empty input files
                        if (isEmpty(f_sample))
                        {
                            continue;
                        }
                    } catch (IOException e)
                    {
                        e.printStackTrace();
                        System.exit(1);
                    }
                    String[] name_split = f_sample.getName().split("\\.");

                    File d_output_sample = extend(d_output.get(), d_group.getName(), d_hm.getName(), name_split[0]);
                    try
                    {
                        makeSureDirectoryExists(d_output_sample);
                    } catch (IOException e)
                    {
                        e.printStackTrace();
                        System.exit(1);
                    }

                    File gtf = f_input_gtf.get();

                    String command_execute = tgeneExecutable.getAbsolutePath();
                    command_execute += " " + f_sample.getAbsolutePath();
                    command_execute += " " + gtf.getAbsolutePath();

                    command_execute += " -oc " + d_output_sample;

                    if (noClosestLocus.get())
                    {
                        command_execute += " --no-closest-locus";
                    }
                    if (noClosestTss.get())
                    {
                        command_execute += " --no-closest-tss";
                    }

                    command_execute += " --max-link-distances " + maxLinkDistance.get();
                    command_execute += " --max-pvalue " + pValue.get();

                    String finalCommand_execute = command_execute;
                    executorService.submit(() -> executeAndWait(finalCommand_execute, logger));
                }
            }
        }
    }
}
