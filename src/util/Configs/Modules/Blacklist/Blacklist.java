package util.Configs.Modules.Blacklist;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;
import java.util.Arrays;
import java.util.List;

public class Blacklist extends AbstractModule
{
    public FileStructure fileStructure;

    public final Config<File> bedFilePath = new Config<>(File.class);
    public final Config<List> signalsToIgnore =
            new Config<>(Arrays.asList("Low_Mappability".toUpperCase(), "High_Signal_Region".toUpperCase()));

    public Blacklist(Config<File> workingDirectory, Config<File> sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
