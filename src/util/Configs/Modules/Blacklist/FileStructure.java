package util.Configs.Modules.Blacklist;

import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.Configs.ConfigTypes.InputFileStructure;
import util.Configs.ConfigTypes.SourceDirectoryFileStructure;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class FileStructure extends AbstractModule
{
    public final GeneratedFileStructure d_root = extend(workingDirectory, "01_blacklisted_regions");
    public final GeneratedFileStructure d_preprocessing = extend(d_root, "01_preprocessing");
    public final GeneratedFileStructure d_preprocessing_perChr = extend(d_preprocessing, "01_perChr");
    public final GeneratedFileStructure d_preprocessing_sorted = extend(d_preprocessing, "02_sorted");
    public final GeneratedFileStructure d_newInput = extend(d_root, "02_new_input");

    public FileStructure(GeneratedFileStructure workingDirectory, SourceDirectoryFileStructure sourceDirectory,
                         Logger logger)
            throws IllegalAccessException, InvocationTargetException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
