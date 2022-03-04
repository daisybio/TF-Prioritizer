package util.Configs.Modules.TpmGcFilterAnalysis;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class FileStructure extends AbstractModule
{
    public final Config<File> d_wd = new Config<>(File.class);

    public final Config<File> d_root = extend(workingDirectory, "A1_TPM_GC_FILTER_ANALYSIS");
    public final Config<File> d_data = extend(d_root, "01_DATA");
    public final Config<File> d_scripts = extend(d_root, "02_RScripts");
    public final Config<File> d_plots = extend(d_root, "03_PLOTS");

    public final Config<String> s_data = new Config<>("data.txt");

    public FileStructure(Config<File> workingDirectory, Config<File> sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
