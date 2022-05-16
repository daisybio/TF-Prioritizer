package lib.ChiPAtlas;

import lib.ExecutableStep;
import org.apache.commons.compress.utils.IOUtils;
import tfprio.tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.TrustAllManager;

import javax.net.ssl.HttpsURLConnection;
import javax.net.ssl.SSLContext;
import javax.net.ssl.TrustManager;
import java.io.*;
import java.net.URL;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import static util.FileManagement.makeSureFileExists;
import static util.ScriptExecution.executeAndWait;

public class GetDataList extends ExecutableStep
{
    private final GeneratedFileStructure f_zipped = TFPRIO.configs.chipAtlas.fileStructure.f_list_zipped;
    private final GeneratedFileStructure f_list = TFPRIO.configs.chipAtlas.fileStructure.f_list_csv;

    private final AbstractConfig<String> genomeVersion = TFPRIO.configs.chipAtlas.genomeVersion;
    private final AbstractConfig<String> tissueType = TFPRIO.configs.chipAtlas.tissueType;
    private final AbstractConfig<String> urlToList = TFPRIO.configs.chipAtlas.urlToList;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>();
    }

    @Override public Set<GeneratedFileStructure> getCreatedFileStructure()
    {
        return new HashSet<>(Arrays.asList(f_zipped, f_list));
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>(Arrays.asList(genomeVersion, tissueType, urlToList));
    }

    @Override protected void execute()
    {
        int maxTries = 10;

        logger.info("Genome version: " + genomeVersion.get());
        logger.info("Tissue type: " + tissueType.get());
        logger.info("Downloading URL: " + urlToList.get());

        makeSureFileExists(f_zipped.get(), logger);

        int i =0;
        boolean worked = false;
        String command = "wget " + urlToList.get() + " -O " + f_zipped.get().getAbsolutePath();
        logger.debug("Executing command: " + command);

        while(i < maxTries && !worked)
        {
            try
            {
                Process wget = Runtime.getRuntime().exec(command);
                wget.waitFor();
                worked = true;
            } catch (IOException | InterruptedException e)
            {
                logger.warn("Failed try #" + i+1);
            }
            i++;
            if (!worked)
            {
                try
                {
                    Thread.sleep(5000);
                } catch (InterruptedException e)
                {
                    logger.error(e.getMessage());
                }
            }
        }

        if (!worked)
        {
            logger.error("Could not get data.");
        }

        if (!f_zipped.get().exists())
        {
            logger.error("Something awkward happened.");
        }

        String command_edited = "unzip " + f_zipped.get() + " -d " + f_list.get().getParentFile().getAbsolutePath();
        executeAndWait(command_edited, logger);

        f_zipped.get().delete();
    }
}
