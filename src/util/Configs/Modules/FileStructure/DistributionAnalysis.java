package util.Configs.Modules.FileStructure;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class DistributionAnalysis extends AbstractModule
{
    public final Config<File> root = extend(workingDirectory, "07_DISTRIBUTION_ANALYSIS");
    public final Config<File> analyzedTfs = extend(root, "01_ANALYSED_TFS");
    public final Config<File> tfTgScores = extend(root, "02_TF_TG_SCORES");
    public final Config<File> tfTgScores_backgroundDistribution = extend(tfTgScores, "01_BACKGROUND_DISTRIBUTION");
    public final Config<File> tfTgScores_backgroundDistribution_all =
            extend(tfTgScores_backgroundDistribution, "01_ALL");
    public final Config<File> tfTgScores_backgroundDistribution_hm = extend(tfTgScores_backgroundDistribution, "02_HM");
    public final Config<File> tfTgScores_tfDistribution = extend(tfTgScores, "02_TF_DISTRIBUTION");
    public final Config<File> tfTgScores_tfDistribution_all = extend(tfTgScores_tfDistribution, "01_ALL");
    public final Config<File> tfTgScores_tfDistribution_hm = extend(tfTgScores_tfDistribution, "02_HM");

    public final Config<File> plotsScripts = extend(root, "03_PLOT_SCRIPTS");
    public final Config<File> plotsScripts_all = extend(plotsScripts, "01_ALL");
    public final Config<File> plotsScripts_hm = extend(plotsScripts, "02_HM");

    public final Config<File> plots = extend(root, "04_DISTR_PLOTS");
    public final Config<File> plots_all = extend(plots, "01_ALL");
    public final Config<File> plots_hm = extend(plots, "02_HM");

    public final Config<File> stats = extend(root, "05_STATS");
    public final Config<File> stats_all = extend(stats, "01_ALL");
    public final Config<File> stats_hm = extend(stats, "02_HM");

    public final Config<File> hypergeometricTest = extend(root, "06_HYPERGEOMETRIC_TEST");
    public final Config<File> mwuScripts = extend(root, "07_SCRIPTS_MANN_WHITNEYU_PLOTS");
    public final Config<File> mwuPlots = extend(root, "08_PLOTS_MANN_WHITNEYU_PLOTS");
    public final Config<File> dcg = extend(root, "09_DISCOUNTED_CUMULATIVE_GAIN");

    public final Config<File> logos = extend(root, "10_LOGOS");
    public final Config<File> logos_biophysicalModel = extend(logos, "01_BIOPHYSICAL_MODEL");
    public final Config<File> logos_tfSequence = extend(logos, "02_TF_SEQUENCE");
    public final Config<File> logos_tfBindingSequence = extend(logos, "03_TF_BINDING_SEQUENCE");
    public final Config<File> logos_tfBindingSequence_data = extend(logos, "00_DATA");

    public final Config<File> cooccurrence = extend(root, "12_COOCCURRENCE");


    public DistributionAnalysis(File workingDirectory, File sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
