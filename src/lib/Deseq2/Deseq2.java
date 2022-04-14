package lib.Deseq2;

import lib.ExecutableStep;
import tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.FileFilters.Filters;

import java.io.*;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;

import static util.ScriptExecution.executeAndWait;

public class Deseq2 extends ExecutableStep
{
    private final AbstractConfig<File> d_scripts = TFPRIO.configs.deSeq2.fileStructure.d_rScripts;
    private final AbstractConfig<File> d_output = TFPRIO.configs.deSeq2.fileStructure.d_outputRaw;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(List.of(d_scripts));
    }

    @Override protected Set<AbstractConfig<File>> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(d_output));
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>();
    }

    @Override protected void execute()
    {
        for (File f_script : Objects.requireNonNull(d_scripts.get().listFiles(Filters.fileFilter)))
        {
            executorService.submit(() -> executeAndWait(f_script, logger));
        }
    }
}
