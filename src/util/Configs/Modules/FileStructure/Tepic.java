package util.Configs.Modules.FileStructure;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class Tepic extends AbstractModule
{
    public final Config<File> outputRaw = extend(workingDirectory, "03_A0_TEPIC_output_raw");
    public final Config<File> outputRaw_shuffle = extend(outputRaw, "03_A1_TEPIC_output_raw_shuffle");
    public final Config<File> postprocessing = extend(workingDirectory, "03_B_TEPIC_postprocessing");
    public final Config<File> postprocessing_input = extend(postprocessing, "input");
    public final Config<File> postprocessing_output = extend(postprocessing, "output");
    public final Config<File> postprocessing_tfs = extend(postprocessing, "TFs");
    public final Config<File> postprocessing_openChromatinViolins =
            extend(postprocessing, "open_chromatin_violin_plots");
    public final Config<File> postprocessing_openChromatinViolins_data = extend(postprocessing, "01_data");
    public final Config<File> postprocessing_openChromatinViolins_script = extend(postprocessing, "02_script");
    public final Config<File> postprocessing_openChromatinViolins_plots = extend(postprocessing, "03_plots");

    public Tepic(File workingDirectory, File sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
