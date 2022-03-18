package lib.Deseq2;

import lib.ExecutableStep;
import tfprio.TFPRIO;
import util.Configs.Config;
import util.ExternalScriptException;
import util.FileFilters.Filters;

import java.io.*;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;

import static util.ScriptExecution.executeAndWait;

public class Deseq2 extends ExecutableStep
{
    private final Config<File> d_scripts = TFPRIO.configs.deSeq2.fileStructure.d_rScripts;
    private final Config<File> d_output = TFPRIO.configs.deSeq2.fileStructure.d_outputRaw;

    @Override protected Set<Config<File>> getRequiredFileStructure()
    {
        return new HashSet<>(List.of(d_scripts));
    }

    @Override protected Set<Config<File>> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(d_output));
    }

    @Override protected Set<Config<?>> getRequiredConfigs()
    {
        return new HashSet<>();
    }

    @Override protected void execute()
    {
        for (File f_script : Objects.requireNonNull(d_scripts.get().listFiles(Filters.fileFilter)))
        {
            executorService.submit(() ->
            {
                try
                {
                    executeAndWait(f_script, logger);
                } catch (ExternalScriptException e)
                {
                    logger.error(e.getMessage());
                    System.exit(1);
                }
            });
        }
    }
}
