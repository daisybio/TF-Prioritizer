package util.Configs.Modules;

import util.Configs.Config;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class Tepic extends AbstractModule
{
    public final Config<File> inputDirectory = new Config<>();
    public final Config<File> inputPrevious = new Config<>();
    public final Config<File> inputOriginal = new Config<>();


    public Tepic(File workingDirectory, File sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
