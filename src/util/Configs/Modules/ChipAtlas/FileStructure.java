package util.Configs.Modules.ChipAtlas;

import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.Configs.ConfigTypes.InputFileStructure;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class FileStructure extends AbstractModule
{
    public final GeneratedFileStructure d_root = extend(workingDirectory, "08_CHIP_ATLAS_EVALUATION_PEAK_DATA");

    public final GeneratedFileStructure d_list = extend(d_root, "01_CHIP_ATLAS_LIST");
    public final GeneratedFileStructure f_list_csv = extend(d_list, "chip_atlas_file_list.csv");
    public final GeneratedFileStructure f_list_zipped = extend(d_list, "chip_atlas_file_list.zip");

    public final GeneratedFileStructure d_peakFiles = extend(d_root, "02_CHIP_ATLAS_PEAK_FILES");

    public FileStructure(GeneratedFileStructure workingDirectory, InputFileStructure sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
