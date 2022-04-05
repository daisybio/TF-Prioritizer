package util.Configs.Modules.Jaspar;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class Jaspar extends AbstractModule
{
    public FileStructure fileStructure;

    public final Config<String> downloadUrl =
            new Config<>("https://jaspar.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_redundant_pfms_jaspar.txt",
                    false);

    public final Config<String> apiCall = new Config<>(
            "curl -X GET https://jaspar.genereg.net/api/v1/matrix/{MATRIXID}/ -H " +
                    "'Accept:application/json' -o '{OUTPUTFILE}' --create-dirs", false);
    public final Config<String> logoDownloadUrl =
            new Config<>("https://jaspar.genereg.net/static/logos/all/svg/", false);


    public Jaspar(Config<File> workingDirectory, Config<File> sourceDirectory, Logger logger)
            throws IllegalAccessException, InvocationTargetException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
