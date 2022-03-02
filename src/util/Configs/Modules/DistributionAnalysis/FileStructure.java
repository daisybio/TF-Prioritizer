package util.Configs.Modules.DistributionAnalysis;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class FileStructure extends AbstractModule
{
    public final Config<File> d_root = extend(workingDirectory, "07_DISTRIBUTION_ANALYSIS");

    public final Config<File> d_analyzedTfs = extend(d_root, "01_ANALYSED_TFS");
    public final Config<File> f_analyzedTfs = extend(d_analyzedTfs, "analysable_tfs.csv");

    public final Config<File> d_tfTgScores = extend(d_root, "02_TF_TG_SCORES");
    public final Config<File> d_tfTgScores_backgroundDistribution = extend(d_tfTgScores, "01_BACKGROUND_DISTRIBUTION");
    public final Config<File> d_tfTgScores_backgroundDistribution_all =
            extend(d_tfTgScores_backgroundDistribution, "01_ALL");
    public final Config<File> d_tfTgScores_backgroundDistribution_hm =
            extend(d_tfTgScores_backgroundDistribution, "02_HM");
    public final Config<String> s_tfTgScores_backgroundDistribution_csv = new Config<>("distribution.csv");

    public final Config<File> d_tfTgScores_tfDistribution = extend(d_tfTgScores, "02_TF_DISTRIBUTION");
    public final Config<File> d_tfTgScores_tfDistribution_all = extend(d_tfTgScores_tfDistribution, "01_ALL");
    public final Config<File> d_tfTgScores_tfDistribution_hm = extend(d_tfTgScores_tfDistribution, "02_HM");

    public final Config<File> d_plotsScripts = extend(d_root, "03_PLOT_SCRIPTS");
    public final Config<File> d_plotsScripts_all = extend(d_plotsScripts, "01_ALL");
    public final Config<File> d_plotsScripts_hm = extend(d_plotsScripts, "02_HM");
    public final Config<String> s_plotsScripts_pythonScript = new Config<>("analyze_tf_to_background.py");

    public final Config<File> d_plots = extend(d_root, "04_DISTR_PLOTS");
    public final Config<File> d_plots_all = extend(d_plots, "01_ALL");
    public final Config<File> d_plots_hm = extend(d_plots, "02_HM");

    public final Config<File> d_stats = extend(d_root, "05_STATS");
    public final Config<File> d_stats_all = extend(d_stats, "01_ALL");
    public final Config<File> d_stats_hm = extend(d_stats, "02_HM");
    public final Config<String> s_stats_csv = new Config<>("stats.csv");

    public final Config<File> d_hypergeometricTest = extend(d_root, "06_HYPERGEOMETRIC_TEST");
    public final Config<File> f_hypergeometricTest_script = extend(d_hypergeometricTest, "HYPERGEOMETRIC_TEST.R");
    public final Config<File> f_hypergeometricTest_output = extend(d_hypergeometricTest, "STATS.csv");

    public final Config<File> d_mwuScripts = extend(d_root, "07_SCRIPTS_MANN_WHITNEYU_PLOTS");
    public final Config<String> s_mwuScripts_script = new Config<>("mann_whitneyU_plots.R");

    public final Config<File> d_mwuPlots = extend(d_root, "08_PLOTS_MANN_WHITNEYU_PLOTS");

    public final Config<File> d_dcg = extend(d_root, "09_DISCOUNTED_CUMULATIVE_GAIN");
    public final Config<File> f_dcg_stats = extend(d_dcg, "dcg_stats.csv");


    public final Config<File> d_logos = extend(d_root, "10_LOGOS");

    public final Config<File> d_logos_biophysicalModel = extend(d_logos, "01_BIOPHYSICAL_MODEL");
    public final Config<String> s_logos_biophysicalModel_data = new Config<>("_energy_matrix.csv");
    public final Config<String> s_logos_biophysicalModel_script = new Config<>("_create_biophysical_model.py");
    public final Config<String> s_logos_biophysicalModel_image = new Config<>("_biophysical_model.png");

    public final Config<File> d_logos_tfSequence = extend(d_logos, "02_TF_SEQUENCE");
    public final Config<File> d_logos_tfSequence_jaspar = extend(d_logos_tfSequence, "00_JASPAR");
    public final Config<File> f_logos_tfSequence_jaspar_pfms = extend(d_logos_tfSequence, "jaspar_pfms.txt");
    public final Config<String> s_logos_tfSequence_jaspar_bash = new Config<>("_curl.sh");
    public final Config<String> s_logos_tfSequence_jaspar_json = new Config<>(".json");
    public final Config<String> s_logos_tfSequence_jaspar_image = new Config<>(".svg");

    public final Config<File> d_logos_tfBindingSequence = extend(d_logos, "03_TF_BINDING_SEQUENCE");
    public final Config<String> s_logos_tfBindingSequence_fasta = new Config<>(".fa");
    public final Config<String> s_logos_tfBindingSequence_script = new Config<>("calc_frequencies.R");
    public final Config<String> s_logos_tfBindingSequence_motif = new Config<>("_frequencies.motif");

    public final Config<File> d_logos_tfBindingSequence_data = extend(d_logos, "00_DATA");

    public final Config<File> d_heatmaps = extend(d_root, "11_DCG_TARGET_GENES_HEATMAP");
    public final Config<File> f_heatmaps_script = extend(d_heatmaps, "heatmap.R");

    public final Config<File> d_cooccurrence = extend(d_root, "12_COOCCURRENCE");
    public final Config<File> f_cooccurrence_mergedBed = extend(d_cooccurrence, "cooccurence_merged.bed");
    public final Config<File> f_cooccurrence_mergedBedSorted = extend(d_cooccurrence, "cooccurence_merged_sorted.bed");
    public final Config<File> f_cooccurrence_script = extend(d_cooccurrence, "sort_and_merge.sh");
    public final Config<File> f_cooccurrence_frequencies = extend(d_cooccurrence, "cooccurence_frequencies.csv");


    public FileStructure(Config<File> workingDirectory, Config<File> sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
