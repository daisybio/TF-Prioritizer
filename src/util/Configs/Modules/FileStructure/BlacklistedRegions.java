package util.Configs.Modules.FileStructure;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class BlacklistedRegions extends AbstractModule
{
    public final Config<File> root =
            new Config<>(new File(workingDirectory.getAbsolutePath() + File.separator + "01_blacklisted_regions"));
    public final Config<File> preprocessing = extend(root, "01_blacklisted_regions");
    public final Config<File> preprocessing_perChr = extend(preprocessing, "01_perChr");
    public final Config<File> preprocessing_sorted = extend(preprocessing, "02_sorted");
    public final Config<File> newInput = extend(root, "02_new_input");

    public BlacklistedRegions(File workingDirectory, File sourceDirectory, Logger logger)
            throws IllegalAccessException, InvocationTargetException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
