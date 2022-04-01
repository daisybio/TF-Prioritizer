package util.Configs.Modules.DistributionAnalysis;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class DistributionAnalysis extends AbstractModule
{
    public FileStructure fileStructure;

    public final Config<Boolean> performAllAnalysis = new Config<>(false);
    public final Config<String> allName = new Config<>("ALL");

    public DistributionAnalysis(Config<File> workingDirectory, Config<File> sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
