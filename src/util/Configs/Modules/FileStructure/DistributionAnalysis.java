package util.Configs.Modules.FileStructure;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class DistributionAnalysis extends AbstractModule
{
    public final Config<File> d_root = extend(workingDirectory, "07_DISTRIBUTION_ANALYSIS");
    public final Config<File> d_analyzedTfs = extend(d_root, "01_ANALYSED_TFS");
    public final Config<File> d_tfTgScores = extend(d_root, "02_TF_TG_SCORES");
    public final Config<File> d_tfTgScores_backgroundDistribution = extend(d_tfTgScores, "01_BACKGROUND_DISTRIBUTION");
    public final Config<File> d_tfTgScores_backgroundDistribution_all =
            extend(d_tfTgScores_backgroundDistribution, "01_ALL");
    public final Config<File> d_tfTgScores_backgroundDistribution_hm =
            extend(d_tfTgScores_backgroundDistribution, "02_HM");
    public final Config<File> d_tfTgScores_tfDistribution = extend(d_tfTgScores, "02_TF_DISTRIBUTION");
    public final Config<File> d_tfTgScores_tfDistribution_all = extend(d_tfTgScores_tfDistribution, "01_ALL");
    public final Config<File> d_tfTgScores_tfDistribution_hm = extend(d_tfTgScores_tfDistribution, "02_HM");

    public final Config<File> d_plotsScripts = extend(d_root, "03_PLOT_SCRIPTS");
    public final Config<File> d_plotsScripts_all = extend(d_plotsScripts, "01_ALL");
    public final Config<File> d_plotsScripts_hm = extend(d_plotsScripts, "02_HM");

    public final Config<File> d_plots = extend(d_root, "04_DISTR_PLOTS");
    public final Config<File> d_plots_all = extend(d_plots, "01_ALL");
    public final Config<File> d_plots_hm = extend(d_plots, "02_HM");

    public final Config<File> d_stats = extend(d_root, "05_STATS");
    public final Config<File> d_stats_all = extend(d_stats, "01_ALL");
    public final Config<File> d_stats_hm = extend(d_stats, "02_HM");

    public final Config<File> d_hypergeometricTest = extend(d_root, "06_HYPERGEOMETRIC_TEST");
    public final Config<File> d_mwuScripts = extend(d_root, "07_SCRIPTS_MANN_WHITNEYU_PLOTS");
    public final Config<File> d_mwuPlots = extend(d_root, "08_PLOTS_MANN_WHITNEYU_PLOTS");
    public final Config<File> d_dcg = extend(d_root, "09_DISCOUNTED_CUMULATIVE_GAIN");

    public final Config<File> d_logos = extend(d_root, "10_LOGOS");
    public final Config<File> d_logos_biophysicalModel = extend(d_logos, "01_BIOPHYSICAL_MODEL");
    public final Config<File> d_logos_tfSequence = extend(d_logos, "02_TF_SEQUENCE");
    public final Config<File> d_logos_tfBindingSequence = extend(d_logos, "03_TF_BINDING_SEQUENCE");
    public final Config<File> d_logos_tfBindingSequence_data = extend(d_logos, "00_DATA");

    public final Config<File> d_cooccurrence = extend(d_root, "12_COOCCURRENCE");


    public DistributionAnalysis(File workingDirectory, File sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
