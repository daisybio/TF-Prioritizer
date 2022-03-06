package util.Configs.Modules.Report;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class OutputStructure extends AbstractModule
{
    public final Config<File> root = extend(workingDirectory, "REPORT");
    public final Config<File> f_home = extend(root, "HOME.html");
    public final Config<File> f_overview = extend(root, "OVERVIEW.html");
    public final Config<File> f_parameters = extend(root, "PARAMETERS.html");
    public final Config<File> f_style = extend(root, "style.css");
    public final Config<File> f_script = extend(root, "script.js");

    public final Config<File> d_validation = extend(root, "VALIDATION");
    public final Config<File> d_validation_logos = extend(d_validation, "LOGOS");
    public final Config<String> s_validation_logos_tfSequence = new Config<>("logosTfSequence");
    public final Config<String> s_validation_logos_tfBindingSequence = new Config<>("logosTfBindingSequence");
    public final Config<String> s_validation_logos_biophysicalModel = new Config<>("logosBiophysicalModel");
    public final Config<String> s_validation_logos_biophysicalModel_png = new Config<>("biophysicalModel.png");


    public final Config<File> d_distribution = extend(root, "DISTRIBUTION");

    public final Config<File> d_regression = extend(root, "REGRESSION");
    public final Config<File> d_regression_performance = extend(d_regression, "Performance");
    public final Config<File> d_regression_performance_barplots = extend(d_regression_performance, "barplots");
    public final Config<File> d_regression_performance_foldChanges = extend(d_regression_performance, "foldChanges");
    public final Config<File> d_regression_performance_heatmap = extend(d_regression_performance, "heatmap");
    public final Config<File> f_regression_performance_html =
            extend(d_regression_performance, "PERFORMANCE_ANALYSIS.html");

    public final Config<File> d_importantLoci = extend(root, "importantLoci");
    public final Config<File> f_importantLoci_html = extend(root, "IMPORTANT_LOCI.html");

    public final Config<File> d_topLog2fc = extend(root, "topLog2fc");
    public final Config<File> f_topLog2fc_html = extend(root, "TOP_LOG2FC.html");

    public final Config<File> d_media = extend(root, "MEDIA");

    public final Config<File> f_cooccurrence_html = extend(root, "COOCCURRENCE.html");

    public OutputStructure(Config<File> workingDirectory, Config<File> sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
