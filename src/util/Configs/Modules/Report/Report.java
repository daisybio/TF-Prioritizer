package util.Configs.Modules.Report;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class Report extends AbstractModule
{
    public FileStructure fileStructure;

    public final Config<String> genecardsUrl =
            new Config<>("https://www.genecards.org/cgi-bin/carddisp" + ".pl?gene={GENE}");


    public Report(Config<File> workingDirectory, Config<File> sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
