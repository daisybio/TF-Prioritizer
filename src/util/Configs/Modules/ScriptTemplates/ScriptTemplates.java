package util.Configs.Modules.ScriptTemplates;

import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.Configs.ConfigTypes.SourceDirectoryFileStructure;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class ScriptTemplates extends AbstractModule
{
    public final SourceDirectoryFileStructure d_root = extend(sourceDirectory, "scripts");
    public final SourceDirectoryFileStructure f_heatmaps = extend(d_root, "heatmap.R");
    public final SourceDirectoryFileStructure f_mapping = extend(d_root, "mapping.R");
    public final SourceDirectoryFileStructure f_deseq2 = extend(d_root, "deseq2.R");
    public final SourceDirectoryFileStructure f_deseq2PreprocessingSingle =
            extend(d_root, "deseq2_preprocessing_single.py");
    public final SourceDirectoryFileStructure f_deseq2PreprocessingCombined =
            extend(d_root, "deseq2_preprocessing_combined.py");
    public final SourceDirectoryFileStructure f_deseq2PreprocessingGetGeneLengths = extend(d_root, "getGeneLengths.R");
    public final SourceDirectoryFileStructure f_deseq2PreprocessingTpm = extend(d_root, "deseq2_tpm.py");
    public final SourceDirectoryFileStructure f_deseq2PreprocessingGetGenePositions =
            extend(d_root, "getGenePositions.R");
    public final SourceDirectoryFileStructure f_deseq2PreprocessingUplift = extend(d_root, "deseq2_upliftPositions.py");
    public final SourceDirectoryFileStructure f_deseq2PreprocessingDbUplift =
            extend(d_root, "deseq2_dbUpliftPositions.py");
    public final SourceDirectoryFileStructure f_plots_openChromatinViolinPlots =
            extend(d_root, "openChromatinViolinPlots.R");
    public final SourceDirectoryFileStructure f_plots_groupPlots = extend(d_root, "groupPlots.py");
    public final SourceDirectoryFileStructure f_distributionPlots = extend(d_root, "distributionPlots.py");
    public final SourceDirectoryFileStructure f_distributionMwuPlots = extend(d_root, "distributionMwuPlots.R");

    public final SourceDirectoryFileStructure f_logos_biophysicalModel = extend(d_root, "logos_biophysicalModel.py");
    public final SourceDirectoryFileStructure f_logos_tfBindingSequence = extend(d_root, "logos_bindingSequence.R");

    public final SourceDirectoryFileStructure f_logDockerStats = extend(d_root, "logDockerStats.sh");

    public final SourceDirectoryFileStructure f_metrics = extend(d_root, "metrics.R");

    public ScriptTemplates(GeneratedFileStructure workingDirectory, SourceDirectoryFileStructure sourceDirectory,
                           Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
