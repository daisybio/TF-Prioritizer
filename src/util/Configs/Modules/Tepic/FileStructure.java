package util.Configs.Modules.Tepic;

import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.Configs.ConfigTypes.InternalConfig;
import util.Configs.ConfigTypes.SourceDirectoryFileStructure;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class FileStructure extends AbstractModule
{
    public final GeneratedFileStructure d_outputRaw = extend(workingDirectory, "03_A0_TEPIC_output_raw");
    public final InternalConfig<String> s_outputRaw_regionsToTargetGenes =
            new InternalConfig<>("regions_to_target_genes.csv");
    public final InternalConfig<String> s_outputRaw_trapSequences = new InternalConfig<>("trap_sequences.csv");

    public final GeneratedFileStructure d_outputRaw_shuffle = extend(d_outputRaw, "03_A1_TEPIC_output_raw_shuffle");
    public final GeneratedFileStructure d_postprocessing = extend(workingDirectory, "03_B_TEPIC_postprocessing");
    public final GeneratedFileStructure f_allTfs = extend(d_postprocessing, "ALL_TFs.csv");

    public final GeneratedFileStructure d_postprocessing_input = extend(d_postprocessing, "input");
    public final GeneratedFileStructure d_postprocessing_output = extend(d_postprocessing, "output");

    public final GeneratedFileStructure d_postprocessing_trapPredictedBeds =
            extend(d_postprocessing, "trap_predicted_tf_bindings");

    public final InternalConfig<String> s_postprocessing_output_meanAffinitiesDir =
            new InternalConfig<>("Mean_Affinities");
    public final InternalConfig<String> s_postprocessing_output_ratiosDir = new InternalConfig<>("Ratio");
    public final InternalConfig<String> s_postprocessing_output_meanAffinitiesFile =
            new InternalConfig<>("Mean_Affinities_");
    public final InternalConfig<String> s_postprocessing_output_ratiosFile = new InternalConfig<>("Ratio_Affinities_");
    public final GeneratedFileStructure d_postprocessing_tfs = extend(d_postprocessing, "TFs");
    public final GeneratedFileStructure f_postprocessing_tfs_csv = extend(d_postprocessing_tfs, "tfs.csv");

    public final GeneratedFileStructure d_postprocessing_openChromatinViolins =
            extend(d_postprocessing, "open_chromatin_violin_plots");
    public final GeneratedFileStructure d_postprocessing_openChromatinViolins_data =
            extend(d_postprocessing_openChromatinViolins, "01_data");
    public final GeneratedFileStructure f_postprocessing_openChromatinViolins_data_csv =
            extend(d_postprocessing_openChromatinViolins_data, "hm_to_open_chromatin_lengths.csv");

    public final GeneratedFileStructure d_postprocessing_openChromatinViolins_script =
            extend(d_postprocessing_openChromatinViolins, "02_script");
    public final GeneratedFileStructure f_postprocessing_openChromatinViolins_script_R =
            extend(d_postprocessing_openChromatinViolins_script, "violin_plots_hm_to_open_chromatin_lengths.R");

    public final GeneratedFileStructure d_postprocessing_openChromatinViolins_plots =
            extend(d_postprocessing_openChromatinViolins, "03_plots");
    public final GeneratedFileStructure f_postprocessing_openChromatinViolins_plots_image =
            extend(d_postprocessing_openChromatinViolins_plots, "open_chromatin_lengths_violin.png");

    public FileStructure(GeneratedFileStructure workingDirectory, SourceDirectoryFileStructure sourceDirectory,
                         Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
