package util.Configs.Modules.Report;

import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.Configs.ConfigTypes.InputFileStructure;
import util.Configs.ConfigTypes.InternalConfig;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class Report extends AbstractModule
{
    public OutputStructure outputStructure;
    public InputStructure inputStructure;

    public final InternalConfig<String> genecardsUrl =
            new InternalConfig<>("https://www.genecards.org/cgi-bin/carddisp" + ".pl?gene={GENE}");


    public Report(GeneratedFileStructure workingDirectory, InputFileStructure sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
