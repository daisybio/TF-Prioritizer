package util.Configs.Modules.Jaspar;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class FileStructure extends AbstractModule
{

    public FileStructure(Config<File> workingDirectory, Config<File> sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
