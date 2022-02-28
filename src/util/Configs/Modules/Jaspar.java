package util.Configs.Modules;

import util.Configs.Config;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class Jaspar extends AbstractModule
{
    public final Config<String> downloadUrl =
            new Config<>("https://jaspar.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_redundant_pfms_jaspar.txt",
                    false);

    public final Config<String> apiCallStart = new Config<>("https://jaspar.genereg.net/api/v1/matrix/", false);
    public final Config<String> apiCallEnd = new Config<>(" -H 'Accept: application/json'", false);
    public final Config<String> logoDownloadUrl =
            new Config<>("https://jaspar.genereg.net/static/logos/all/svg/", false);


    public Jaspar(File workingDirectory, File sourceDirectory, Logger logger)
            throws IllegalAccessException, InvocationTargetException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
