package lib.AngularReport;

import lib.ExecutableStep;
import tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import static util.FileManagement.*;

public class PostProcessing extends ExecutableStep
{
    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>();
    }

    @Override public Set<GeneratedFileStructure> getCreatedFileStructure()
    {
        return new HashSet<>();
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>();
    }

    @Override protected void execute()
    {
        File f_index = extend(TFPRIO.configs.angularReport.fileStructure.d_output.get(), "index.html");
        try
        {
            String content = readFile(f_index);
            writeFile(f_index, content.replace(" type=\"module\"", ""));
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }
    }
}
