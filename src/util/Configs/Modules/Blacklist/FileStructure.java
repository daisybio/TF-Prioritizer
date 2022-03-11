package util.Configs.Modules.Blacklist;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class FileStructure extends AbstractModule
{
    public final Config<File> d_root = extend(workingDirectory, "01_blacklisted_regions");
    public final Config<File> d_preprocessing = extend(d_root, "01_preprocessing");
    public final Config<File> d_preprocessing_perChr = extend(d_preprocessing, "01_perChr");
    public final Config<File> d_preprocessing_sorted = extend(d_preprocessing, "02_sorted");
    public final Config<File> d_newInput = extend(d_root, "02_new_input");

    public FileStructure(Config<File> workingDirectory, Config<File> sourceDirectory, Logger logger)
            throws IllegalAccessException, InvocationTargetException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
