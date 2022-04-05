package util.Configs.Modules.ScriptTemplates;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class ScriptTemplates extends AbstractModule
{
    public final Config<File> d_root = extend(sourceDirectory, "scripts");
    public final Config<File> f_heatmaps = extend(d_root, "heatmap.R");
    public final Config<File> f_mapping = extend(d_root, "mapping.R");
    public final Config<File> f_deseq2 = extend(d_root, "deseq2.R");
    public final Config<File> f_deseq2PreprocessingSingle = extend(d_root, "deseq2_preprocessing_single.py");
    public final Config<File> f_deseq2PreprocessingCombined = extend(d_root, "deseq2_preprocessing_combined.py");
    public final Config<File> f_deseq2PreprocessingGetGeneLengths = extend(d_root, "getGeneLengths.R");
    public final Config<File> f_deseq2PreprocessingTpm = extend(d_root, "deseq2_tpm.py");
    public final Config<File> f_deseq2PreprocessingGetGenePositions = extend(d_root, "getGenePositions.R");
    public final Config<File> f_deseq2PreprocessingUplift = extend(d_root, "deseq2_upliftPositions.py");
    public final Config<File> f_plots_openChromatinViolinPlots = extend(d_root, "openChromatinViolinPlots.R");
    public final Config<File> f_plots_groupPlots = extend(d_root, "groupPlots.py");
    public final Config<File> f_distributionPlots = extend(d_root, "distributionPlots.py");
    public final Config<File> f_distributionMwuPlots = extend(d_root, "distributionMwuPlots.R");

    public final Config<File> f_logos_biophysicalModel = extend(d_root, "logos_biophysicalModel.py");

    public ScriptTemplates(Config<File> workingDirectory, Config<File> sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
