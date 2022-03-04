package util.Configs.Modules.Tepic;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class FileStructure extends AbstractModule
{
    public final Config<File> d_outputRaw = extend(workingDirectory, "03_A0_TEPIC_output_raw");
    public final Config<String> s_outputRaw_regionsToTargetGenes = new Config<>("regions_to_target_genes.csv");
    public final Config<String> s_outputRaw_trapSequences = new Config<>("trap_sequences.csv");

    public final Config<File> d_outputRaw_shuffle = extend(d_outputRaw, "03_A1_TEPIC_output_raw_shuffle");
    public final Config<File> d_postprocessing = extend(workingDirectory, "03_B_TEPIC_postprocessing");
    public final Config<File> f_allTfs = extend(d_postprocessing, "ALL_TFs.csv");

    public final Config<File> d_postprocessing_input = extend(d_postprocessing, "input");
    public final Config<File> d_postprocessing_output = extend(d_postprocessing, "output");
    public final Config<String> s_postprocessing_output_meanAffinitiesDir = new Config<>("Mean_Affinities");
    public final Config<String> s_postprocessing_output_ratiosDir = new Config<>("Ratio");
    public final Config<String> s_postprocessing_output_meanAffinitiesFile = new Config<>("Mean_Affinities_");
    public final Config<String> s_postprocessing_output_ratiosFile = new Config<>("Ratio_Affinities_");
    public final Config<File> d_postprocessing_tfs = extend(d_postprocessing, "TFs");
    public final Config<File> f_postprocessing_tfs_csv = extend(d_postprocessing_tfs, "tfs.csv");

    public final Config<File> d_postprocessing_openChromatinViolins =
            extend(d_postprocessing, "open_chromatin_violin_plots");
    public final Config<File> d_postprocessing_openChromatinViolins_data = extend(d_postprocessing, "01_data");
    public final Config<File> f_postprocessing_openChromatinViolins_data_csv =
            extend(d_postprocessing_openChromatinViolins_data, "hm_to_open_chromatin_lengths.csv");

    public final Config<File> d_postprocessing_openChromatinViolins_script = extend(d_postprocessing, "02_script");
    public final Config<File> f_postprocessing_openChromatinViolins_script_R =
            extend(d_postprocessing_openChromatinViolins_script, "violin_plots_hm_to_open_chromatin_lengths.R");

    public final Config<File> d_postprocessing_openChromatinViolins_plots = extend(d_postprocessing, "03_plots");
    public final Config<File> f_postprocessing_openChromatinViolins_plots_image =
            extend(d_postprocessing_openChromatinViolins_plots, "open_chromatin_lengths_violin.png");

    public FileStructure(Config<File> workingDirectory, Config<File> sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
