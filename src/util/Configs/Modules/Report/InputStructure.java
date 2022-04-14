package util.Configs.Modules.Report;

import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.Configs.ConfigTypes.InputFileStructure;
import util.Configs.ConfigTypes.SourceDirectoryFileStructure;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class InputStructure extends AbstractModule
{
    public final SourceDirectoryFileStructure d_root = extend(sourceDirectory, "ext" + File.separator + "REPORT");
    public final SourceDirectoryFileStructure d_media = extend(d_root, "MEDIA");
    public final SourceDirectoryFileStructure f_importantLoci = extend(d_root, "IMPORTANT_LOCI.html");
    public final SourceDirectoryFileStructure f_topLog2fc = extend(d_root, "TOP_LOG2FC.html");
    public final SourceDirectoryFileStructure f_cooccurrence = extend(d_root, "COOCCURRENCE.html");
    public final SourceDirectoryFileStructure f_overview = extend(d_root, "OVERVIEW.html");
    public final SourceDirectoryFileStructure f_frame = extend(d_root, "FRAME.html");
    public final SourceDirectoryFileStructure f_style = extend(d_root, "style.css");
    public final SourceDirectoryFileStructure f_script = extend(d_root, "script.js");

    public final SourceDirectoryFileStructure d_basicData = extend(d_root, "BASICDATA");
    public final SourceDirectoryFileStructure f_basicData = extend(d_basicData, "BASICDATA.html");
    public final SourceDirectoryFileStructure f_basicDataEntry = extend(d_basicData, "ENTRY.html");
    public final SourceDirectoryFileStructure f_basicDataGeneID = extend(d_basicData, "GENEID.html");

    public final SourceDirectoryFileStructure d_home = extend(d_root, "HOME");
    public final SourceDirectoryFileStructure f_home_buttonbar = extend(d_home, "BUTTONBAR.html");
    public final SourceDirectoryFileStructure f_home = extend(d_home, "HOME.html");
    public final SourceDirectoryFileStructure f_home_tf = extend(d_home, "TF.html");
    public final SourceDirectoryFileStructure f_home_tfGroup = extend(d_home, "TF_GROUP.html");

    public final SourceDirectoryFileStructure d_parameters = extend(d_root, "PARAMETERS");
    public final SourceDirectoryFileStructure f_parameters = extend(d_parameters, "PARAMETERS.html");
    public final SourceDirectoryFileStructure f_parameters_tool = extend(d_parameters, "TOOL.html");
    public final SourceDirectoryFileStructure f_parameters_parameter = extend(d_parameters, "PARAMETER.html");

    public final SourceDirectoryFileStructure d_validation = extend(d_root, "VALIDATION");
    public final SourceDirectoryFileStructure f_validation = extend(d_validation, "VALIDATION.html");

    public final SourceDirectoryFileStructure d_distribution = extend(d_root, "DISTRIBUTION");
    public final SourceDirectoryFileStructure f_distribution = extend(d_distribution, "DISTRIBUTION.html");

    public final SourceDirectoryFileStructure d_regression = extend(d_root, "REGRESSION");
    public final SourceDirectoryFileStructure f_regression = extend(d_regression, "REGRESSION.html");
    public final SourceDirectoryFileStructure f_regressionPerformance =
            extend(d_regression, "PERFORMANCE_ANALYSIS.html");

    public InputStructure(GeneratedFileStructure workingDirectory, SourceDirectoryFileStructure sourceDirectory,
                          Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
