package util.Configs.Modules.Jaspar;

import util.Configs.ConfigTypes.*;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class Jaspar extends AbstractModule
{
    public FileStructure fileStructure;

    public final InternalConfig<String> downloadUrl = new InternalConfig<>(
            "https://jaspar.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_redundant_pfms_jaspar.txt");

    public final InternalConfig<String> apiCall = new InternalConfig<>(
            "curl -X GET https://jaspar.genereg.net/api/v1/matrix/{MATRIXID}/ -H " +
                    "'Accept:application/json' -o '{OUTPUTFILE}'");
    public final InternalConfig<String> logoDownloadUrl =
            new InternalConfig<>("https://jaspar.genereg.net/static/logos/all/svg/");


    public Jaspar(GeneratedFileStructure workingDirectory, SourceDirectoryFileStructure sourceDirectory, Logger logger)
            throws IllegalAccessException, InvocationTargetException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
