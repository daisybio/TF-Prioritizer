package util.Configs.Modules.Report;

import util.Configs.ConfigTypes.*;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class OutputStructure extends AbstractModule
{
    public final GeneratedFileStructure root = extend(workingDirectory, "REPORT");
    public final GeneratedFileStructure f_home = extend(root, "HOME.html");
    public final GeneratedFileStructure f_overview = extend(root, "OVERVIEW.html");
    public final GeneratedFileStructure f_parameters = extend(root, "PARAMETERS.html");
    public final GeneratedFileStructure f_style = extend(root, "style.css");
    public final GeneratedFileStructure f_script = extend(root, "script.js");

    public final GeneratedFileStructure d_validation = extend(root, "VALIDATION");
    public final GeneratedFileStructure d_validation_logos = extend(d_validation, "LOGOS");
    public final InternalConfig<String> s_validation_logos_tfSequence = new InternalConfig<>("logosTfSequence");
    public final InternalConfig<String> s_validation_logos_tfBindingSequence =
            new InternalConfig<>("logosTfBindingSequence");
    public final InternalConfig<String> s_validation_logos_biophysicalModel =
            new InternalConfig<>("logosBiophysicalModel");
    public final InternalConfig<String> s_validation_logos_biophysicalModel_png =
            new InternalConfig<>("biophysicalModel.png");


    public final GeneratedFileStructure d_distribution = extend(root, "DISTRIBUTION");

    public final GeneratedFileStructure d_regression = extend(root, "REGRESSION");
    public final GeneratedFileStructure d_regression_performance = extend(d_regression, "Performance");
    public final GeneratedFileStructure d_regression_performance_barplots =
            extend(d_regression_performance, "barplots");
    public final GeneratedFileStructure d_regression_performance_foldChanges =
            extend(d_regression_performance, "foldChanges");
    public final GeneratedFileStructure d_regression_performance_heatmap = extend(d_regression_performance, "heatmap");
    public final GeneratedFileStructure f_regression_performance_html =
            extend(d_regression_performance, "PERFORMANCE_ANALYSIS.html");

    public final GeneratedFileStructure d_importantLoci = extend(root, "importantLoci");
    public final GeneratedFileStructure f_importantLoci_html = extend(root, "IMPORTANT_LOCI.html");

    public final GeneratedFileStructure d_topLog2fc = extend(root, "topLog2fc");
    public final GeneratedFileStructure f_topLog2fc_html = extend(root, "TOP_LOG2FC.html");

    public final GeneratedFileStructure d_media = extend(root, "MEDIA");

    public final GeneratedFileStructure f_cooccurrence_html = extend(root, "COOCCURRENCE.html");

    public OutputStructure(GeneratedFileStructure workingDirectory, InputFileStructure sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
