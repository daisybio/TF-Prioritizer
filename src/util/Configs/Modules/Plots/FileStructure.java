package util.Configs.Modules.Plots;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class FileStructure extends AbstractModule
{
    public final Config<File> d_data = extend(workingDirectory, "06_A_PLOTS_DATA");
    public final Config<String> s_data_hmLevelDifferent = new Config<>("all_data_different.csv");
    public final Config<String> s_data_hmLevelSame = new Config<>("all_data_same.csv");

    public final Config<File> d_output = extend(workingDirectory, "06_B_PLOTS_output");

    public final Config<File> d_analysisData = extend(workingDirectory, "06_C_PLOTS_ANALYSIS");
    public final Config<File> d_analysisData_tpLevel = extend(d_analysisData, "01_TP_LEVEL");
    public final Config<File> d_analysisData_hmLevel = extend(d_analysisData, "02_HM_LEVEL");
    public final Config<File> d_analysisData_websiteOverview = extend(d_analysisData, "03_WEBSITE_OVERVIEW");
    public final Config<String> s_analysisData_websiteOverview_tfAvailable = new Config<>("available_tfs.csv");

    public final Config<File> d_analysisData_websiteNumberTfs = extend(d_analysisData, "04_WEBSITE_NUMBER_TFs");

    public final Config<File> d_targetGenes = extend(workingDirectory, "06_D_PLOTS_TARGET_GENES");
    public final Config<File> d_targetGenes_allDifferent = extend(d_targetGenes, "all_data_different");
    public final Config<File> d_targetGenes_allSame = extend(d_targetGenes, "all_data_same");

    public final Config<File> d_targetGenesDcg = extend(workingDirectory, "06_E_PLOTS_TARGET_GENES_DCG");

    public FileStructure(Config<File> workingDirectory, Config<File> sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
