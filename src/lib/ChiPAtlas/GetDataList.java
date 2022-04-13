package lib.ChiPAtlas;

import lib.ExecutableStep;
import org.apache.commons.compress.utils.IOUtils;
import tfprio.TFPRIO;
import util.Configs.Config;
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
    private final Config<File> f_zipped = TFPRIO.configs.chipAtlas.fileStructure.f_list_zipped;
    private final Config<File> f_list = TFPRIO.configs.chipAtlas.fileStructure.f_list_csv;

    private final Config<String> genomeVersion = TFPRIO.configs.chipAtlas.genomeVersion;
    private final Config<String> tissueType = TFPRIO.configs.chipAtlas.tissueType;
    private final Config<String> urlToList = TFPRIO.configs.chipAtlas.urlToList;

    @Override protected Set<Config<File>> getRequiredFileStructure()
    {
        return new HashSet<>();
    }

    @Override protected Set<Config<File>> getCreatedFileStructure()
    {
        return new HashSet<>(Arrays.asList(f_zipped, f_list));
    }

    @Override protected Set<Config<?>> getRequiredConfigs()
    {
        return new HashSet<>(Arrays.asList(genomeVersion, tissueType, urlToList));
    }

    @Override protected void execute()
    {
        logger.info("Genome version: " + genomeVersion.get());
        logger.info("Tissue type: " + tissueType.get());
        logger.info("Downloading URL: " + urlToList.get());

        TrustManager[] trustAllCerts = new TrustManager[]{new TrustAllManager()};
        try
        {
            SSLContext sc = SSLContext.getInstance("SSL");
            sc.init(null, trustAllCerts, new java.security.SecureRandom());
            HttpsURLConnection.setDefaultSSLSocketFactory(sc.getSocketFactory());
        } catch (Exception e)
        {
            logger.error(e.getMessage());
        }

        makeSureFileExists(f_zipped.get(), logger);
        try (InputStream inputStream = new URL(urlToList.get()).openConnection().getInputStream();
             OutputStream outputStream = new FileOutputStream(f_zipped.get()))
        {
            IOUtils.copy(inputStream, outputStream);
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        String command_edited = "unzip " + f_zipped.get() + " -d " + f_list.get().getParentFile().getAbsolutePath();
        executeAndWait(command_edited, logger);
    }
}
