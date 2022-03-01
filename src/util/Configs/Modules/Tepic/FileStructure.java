package util.Configs.Modules.Tepic;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class FileStructure extends AbstractModule
{
    public final Config<File> d_outputRaw = extend(workingDirectory, "03_A0_TEPIC_output_raw");
    public final Config<File> d_outputRaw_shuffle = extend(d_outputRaw, "03_A1_TEPIC_output_raw_shuffle");
    public final Config<File> d_postprocessing = extend(workingDirectory, "03_B_TEPIC_postprocessing");
    public final Config<File> d_postprocessing_input = extend(d_postprocessing, "input");
    public final Config<File> d_postprocessing_output = extend(d_postprocessing, "output");
    public final Config<File> d_postprocessing_tfs = extend(d_postprocessing, "TFs");
    public final Config<File> d_postprocessing_openChromatinViolins =
            extend(d_postprocessing, "open_chromatin_violin_plots");
    public final Config<File> d_postprocessing_openChromatinViolins_data = extend(d_postprocessing, "01_data");
    public final Config<File> d_postprocessing_openChromatinViolins_script = extend(d_postprocessing, "02_script");
    public final Config<File> d_postprocessing_openChromatinViolins_plots = extend(d_postprocessing, "03_plots");

    public FileStructure(File workingDirectory, File sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
