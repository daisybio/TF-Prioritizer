package util.Configs.Modules.TpmGcFilterAnalysis;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;
import java.util.List;

public class TpmGcFilterAnalysis extends AbstractModule
{
    public FileStructure fileStructure;

    public final Config<List> tfList = new Config<>(List.class);
    public final Config<Boolean> countZeros = new Config<>(true);

    public TpmGcFilterAnalysis(Config<File> workingDirectory, Config<File> sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
