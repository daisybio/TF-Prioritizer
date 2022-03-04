package util.Configs.Modules.Report;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class Input extends AbstractModule
{
    public final Config<File> d_root = extend(sourceDirectory, "ext" + File.separator + "REPORT");
    public final Config<File> d_media = extend(d_root, "MEDIA");
    public final Config<File> f_logo = extend(d_media, "logo.png");
    public final Config<File> f_isAvailable = extend(d_media, "is_available.png");
    public final Config<File> f_notAvailable = extend(d_media, "not_available.png");
    public final Config<File> f_importantLoci = extend(d_root, "IMPORTANT_LOCI.html");
    public final Config<File> f_topLog2fc = extend(d_root, "TOP_LOG2FC.html");
    public final Config<File> f_cooccurrence = extend(d_root, "COOCCURRENCE.html");
    public final Config<File> f_overview = extend(d_root, "OVERVIEW.html");
    public final Config<File> f_frame = extend(d_root, "FRAME.html");
    public final Config<File> f_style = extend(d_root, "style.css");
    public final Config<File> f_script = extend(d_root, "script.js");

    public final Config<File> d_basicData = extend(d_root, "BASICDATA");
    public final Config<File> f_basicData = extend(d_basicData, "BASICDATA.html");
    public final Config<File> f_basicDataEntry = extend(d_basicData, "ENTRY.html");
    public final Config<File> f_basicDataGeneID = extend(d_basicData, "GENEID.html");

    public final Config<File> d_home = extend(d_root, "HOME");
    public final Config<File> f_home_buttonbar = extend(d_home, "BUTTONBAR.html");
    public final Config<File> f_home = extend(d_home, "HOME.html");
    public final Config<File> f_home_tf = extend(d_home, "TF.html");
    public final Config<File> f_home_tfGroup = extend(d_home, "TF_GROUP.html");

    public final Config<File> d_parameters = extend(d_root, "PARAMETERS");
    public final Config<File> f_parameters = extend(d_parameters, "PARAMETERS.html");
    public final Config<File> f_parameters_tool = extend(d_parameters, "TOOL.html");
    public final Config<File> f_parameters_parameter = extend(d_parameters, "PARAMETER.html");

    public final Config<File> d_validation = extend(d_root, "VALIDATION");
    public final Config<File> f_validation = extend(d_validation, "VALIDATION.html");

    public final Config<File> d_distribution = extend(d_root, "DISTRIBUTION");
    public final Config<File> f_distribution = extend(d_distribution, "DISTRIBUTION.html");

    public final Config<File> d_regression = extend(d_root, "REGRESSION");
    public final Config<File> f_regression = extend(d_regression, "REGRESSION.html");
    public final Config<File> f_regressionPerformance = extend(d_regression, "PERFORMANCE_ANALYSIS.html");

    public Input(Config<File> workingDirectory, Config<File> sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
