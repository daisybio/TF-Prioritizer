package util.Configs.Modules.MixOptions;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class MixOptions extends AbstractModule
{
    public FileStructure fileStructure;

    public final Config<String> level = new Config<>(String.class);
    public final Config<String> option = new Config<>(String.class);
    public final Config<Integer> occurrenceIntersection = new Config<>(2);
    public final Config<Boolean> mutuallyExclusive = new Config<>(false);
    public final Config<Boolean> mutuallyExclusiveDifferentialPeakSignals = new Config<>(true);

    public MixOptions(Config<File> workingDirectory, Config<File> sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
