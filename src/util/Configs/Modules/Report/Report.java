package util.Configs.Modules.Report;

import util.Configs.ConfigTypes.*;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.lang.reflect.InvocationTargetException;

public class Report extends AbstractModule
{
    public OutputStructure outputStructure;
    public InputStructure inputStructure;

    public final InternalConfig<String> genecardsUrl =
            new InternalConfig<>("https://www.genecards.org/cgi-bin/carddisp" + ".pl?gene={GENE}");


    public Report(GeneratedFileStructure workingDirectory, SourceDirectoryFileStructure sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
