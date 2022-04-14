package util.Configs.Modules.DistributionAnalysis;

import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.Configs.ConfigTypes.InputFileStructure;
import util.Configs.ConfigTypes.InternalConfig;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class FileStructure extends AbstractModule
{
    public final GeneratedFileStructure d_root = extend(workingDirectory, "07_DISTRIBUTION_ANALYSIS");

    public final GeneratedFileStructure d_analyzedTfs = extend(d_root, "01_ANALYSED_TFS");
    public final GeneratedFileStructure f_analyzedTfs = extend(d_analyzedTfs, "analysable_tfs.csv");

    public final GeneratedFileStructure d_tfTgScores = extend(d_root, "02_TF_TG_SCORES");
    public final GeneratedFileStructure d_tfTgScores_backgroundDistribution =
            extend(d_tfTgScores, "01_BACKGROUND_DISTRIBUTION");
    public final GeneratedFileStructure d_tfTgScores_backgroundDistribution_all =
            extend(d_tfTgScores_backgroundDistribution, "01_ALL");
    public final GeneratedFileStructure d_tfTgScores_backgroundDistribution_hm =
            extend(d_tfTgScores_backgroundDistribution, "02_HM");
    public final InternalConfig<String> s_tfTgScores_backgroundDistribution_csv =
            new InternalConfig<>("distribution.csv");

    public final GeneratedFileStructure d_tfTgScores_tfDistribution = extend(d_tfTgScores, "02_TF_DISTRIBUTION");
    public final GeneratedFileStructure d_tfTgScores_tfDistribution_all = extend(d_tfTgScores_tfDistribution, "01_ALL");
    public final GeneratedFileStructure d_tfTgScores_tfDistribution_hm = extend(d_tfTgScores_tfDistribution, "02_HM");

    public final GeneratedFileStructure d_tfTgScores_raw = extend(d_tfTgScores, "00_RAW");

    public final GeneratedFileStructure d_plotsScripts = extend(d_root, "03_PLOT_SCRIPTS");
    public final GeneratedFileStructure d_plotsScripts_all = extend(d_plotsScripts, "01_ALL");
    public final GeneratedFileStructure d_plotsScripts_hm = extend(d_plotsScripts, "02_HM");
    public final InternalConfig<String> s_plotsScripts_pythonScript =
            new InternalConfig<>("analyze_tf_to_background.py");

    public final GeneratedFileStructure d_plots = extend(d_root, "04_DISTR_PLOTS");
    public final GeneratedFileStructure d_plots_all = extend(d_plots, "01_ALL");
    public final GeneratedFileStructure d_plots_hm = extend(d_plots, "02_HM");

    public final GeneratedFileStructure d_stats = extend(d_root, "05_STATS");
    public final GeneratedFileStructure d_stats_all = extend(d_stats, "01_ALL");
    public final GeneratedFileStructure d_stats_hm = extend(d_stats, "02_HM");
    public final InternalConfig<String> s_stats_csv = new InternalConfig<>("stats.csv");

    public final GeneratedFileStructure d_hypergeometricTest = extend(d_root, "06_HYPERGEOMETRIC_TEST");
    public final GeneratedFileStructure f_hypergeometricTest_script =
            extend(d_hypergeometricTest, "HYPERGEOMETRIC_TEST.R");
    public final GeneratedFileStructure f_hypergeometricTest_output = extend(d_hypergeometricTest, "STATS.csv");

    public final GeneratedFileStructure d_mwuScripts = extend(d_root, "07_SCRIPTS_MANN_WHITNEYU_PLOTS");
    public final GeneratedFileStructure d_mwuScriptsHm = extend(d_mwuScripts, "02_HM");
    public final GeneratedFileStructure d_mwuScriptsAll = extend(d_mwuScripts, "01_ALL");
    public final InternalConfig<String> s_mwuScripts_script = new InternalConfig<>("mann_whitneyU_plots.R");

    public final GeneratedFileStructure d_mwuPlots = extend(d_root, "08_PLOTS_MANN_WHITNEYU_PLOTS");

    public final GeneratedFileStructure d_dcg = extend(d_root, "09_DISCOUNTED_CUMULATIVE_GAIN");
    public final GeneratedFileStructure f_dcg_stats = extend(d_dcg, "dcg_stats.csv");
    public final GeneratedFileStructure d_dcg_targetGenes = extend(d_dcg, "targetGenes");


    public final GeneratedFileStructure d_logos = extend(d_root, "10_LOGOS");

    public final GeneratedFileStructure d_logos_biophysicalModel = extend(d_logos, "01_BIOPHYSICAL_MODEL");
    public final InternalConfig<String> s_logos_biophysicalModel_data = new InternalConfig<>("energy_matrix.csv");
    public final InternalConfig<String> s_logos_biophysicalModel_script = new InternalConfig<>("biophysical_model.py");
    public final InternalConfig<String> s_logos_biophysicalModel_image = new InternalConfig<>("biophysical_model.png");

    public final GeneratedFileStructure d_logos_tfSequence = extend(d_logos, "02_TF_SEQUENCE");
    public final GeneratedFileStructure d_logos_tfSequence_jaspar = extend(d_logos_tfSequence, "00_JASPAR");
    public final GeneratedFileStructure f_logos_tfSequence_jaspar_pfms = extend(d_logos_tfSequence, "jaspar_pfms.txt");
    public final InternalConfig<String> s_logos_tfSequence_jaspar_bash = new InternalConfig<>("_curl.sh");
    public final InternalConfig<String> s_logos_tfSequence_jaspar_json = new InternalConfig<>(".json");
    public final InternalConfig<String> s_logos_tfSequence_jaspar_image = new InternalConfig<>(".svg");

    public final GeneratedFileStructure d_logos_tfBindingSequence = extend(d_logos, "03_TF_BINDING_SEQUENCE");
    public final InternalConfig<String> s_logos_tfBindingSequence_fasta = new InternalConfig<>(".fa");
    public final InternalConfig<String> s_logos_tfBindingSequence_script = new InternalConfig<>("calc_frequencies.R");
    public final InternalConfig<String> s_logos_tfBindingSequence_motif = new InternalConfig<>("_frequencies.motif");

    public final GeneratedFileStructure d_logos_tfBindingSequence_data = extend(d_logos_tfBindingSequence, "00_DATA");

    public final GeneratedFileStructure d_heatmaps = extend(d_root, "11_DCG_TARGET_GENES_HEATMAP");
    public final GeneratedFileStructure f_heatmaps_script = extend(d_heatmaps, "heatmap.R");

    public final GeneratedFileStructure d_cooccurrence = extend(d_root, "12_COOCCURENCE");
    public final GeneratedFileStructure f_cooccurrence_concatenatedBed = extend(d_cooccurrence, "concatenated.bed");
    public final GeneratedFileStructure f_cooccurrence_sortedBed = extend(d_cooccurrence, "sorted.bed");
    public final GeneratedFileStructure f_cooccurrence_mergedBed = extend(d_cooccurrence, "merged.bed");
    public final GeneratedFileStructure f_cooccurrence_frequencies = extend(d_cooccurrence, "frequencies.csv");


    public FileStructure(GeneratedFileStructure workingDirectory, InputFileStructure sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
