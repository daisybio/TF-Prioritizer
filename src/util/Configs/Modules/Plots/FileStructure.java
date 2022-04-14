package util.Configs.Modules.Plots;

import util.Configs.ConfigTypes.*;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class FileStructure extends AbstractModule
{
    public final GeneratedFileStructure d_data = extend(workingDirectory, "06_A_PLOTS_DATA");
    public final InternalConfig<String> s_data_hmLevelDifferent = new InternalConfig<>("all_data_different.csv");
    public final InternalConfig<String> s_data_hmLevelSame = new InternalConfig<>("all_data_same.csv");

    public final GeneratedFileStructure d_output = extend(workingDirectory, "06_B_PLOTS_output");

    public final GeneratedFileStructure d_analysisData = extend(workingDirectory, "06_C_PLOTS_ANALYSIS");
    public final GeneratedFileStructure d_analysisData_tpLevel = extend(d_analysisData, "01_TP_LEVEL");
    public final GeneratedFileStructure d_analysisData_hmLevel = extend(d_analysisData, "02_HM_LEVEL");
    public final GeneratedFileStructure d_analysisData_websiteOverview = extend(d_analysisData, "03_WEBSITE_OVERVIEW");
    public final InternalConfig<String> s_analysisData_websiteOverview_tfAvailable =
            new InternalConfig<>("available_tfs.csv");

    public final GeneratedFileStructure d_analysisData_websiteNumberTfs =
            extend(d_analysisData, "04_WEBSITE_NUMBER_TFs");

    public final GeneratedFileStructure d_targetGenes = extend(workingDirectory, "06_D_PLOTS_TARGET_GENES");
    public final GeneratedFileStructure d_targetGenes_allDifferent = extend(d_targetGenes, "all_data_different");
    public final GeneratedFileStructure d_targetGenes_allSame = extend(d_targetGenes, "all_data_same");

    public final GeneratedFileStructure d_targetGenesDcg = extend(workingDirectory, "06_E_PLOTS_TARGET_GENES_DCG");

    public FileStructure(GeneratedFileStructure workingDirectory, SourceDirectoryFileStructure sourceDirectory,
                         Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
