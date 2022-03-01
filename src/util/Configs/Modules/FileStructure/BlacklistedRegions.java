package util.Configs.Modules.FileStructure;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class BlacklistedRegions extends AbstractModule
{
    public final Config<File> d_root =
            new Config<>(new File(workingDirectory.getAbsolutePath() + File.separator + "01_blacklisted_regions"));
    public final Config<File> d_preprocessing = extend(d_root, "01_blacklisted_regions");
    public final Config<File> d_preprocessing_perChr = extend(d_preprocessing, "01_perChr");
    public final Config<File> d_preprocessing_sorted = extend(d_preprocessing, "02_sorted");
    public final Config<File> d_newInput = extend(d_root, "02_new_input");

    public BlacklistedRegions(File workingDirectory, File sourceDirectory, Logger logger)
            throws IllegalAccessException, InvocationTargetException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
