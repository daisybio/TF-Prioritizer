package util.Configs.Modules.Report;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class Output extends AbstractModule
{
    public final Config<File> root = extend(workingDirectory, "REPORT");
    public final Config<File> f_home = extend(root, "HOME.html");
    public final Config<File> f_overview = extend(root, "OVERVIEW.html");
    public final Config<File> f_parameters = extend(root, "PARAMETERS.html");
    public final Config<File> f_style = extend(root, "style.css");
    public final Config<File> f_script = extend(root, "script.js");

    public final Config<File> d_validation = extend(root, "VALIDATION");
    public final Config<File> d_validation_logos = extend(d_validation, "LOGOS");
    public final Config<File> d_validation_logos_tfSequence = extend(d_validation_logos, "tfSequence");
    public final Config<File> d_validation_logos_tfBindingSequence = extend(d_validation_logos, "tfBindingSequence");
    public final Config<File> d_validation_logos_biophysicalModel = extend(d_validation_logos, "biophysicalModel");
    public final Config<File> d_validation_logos_biophysicalModel_png =
            extend(d_validation_logos_biophysicalModel, "biophysicalModel.png");


    public final Config<File> d_distribution = extend(root, "DISTRIBUTION");

    public final Config<File> d_regression = extend(root, "REGRESSION");

    public final Config<File> d_importantLoci = extend(root, "importantLoci");
    public final Config<File> f_importantLoci_html = extend(d_importantLoci, "IMPORTANT_LOCI.html");

    public final Config<File> d_topLog2fc = extend(root, "topLog2fc");
    public final Config<File> f_topLog2fc_html = extend(d_topLog2fc, "TOP_LOG2FC.html");

    public final Config<File> d_media = extend(root, "MEDIA");

    public final Config<File> f_cooccurrence_html = extend(root, "COOCCURRENCE.html");

    public Output(File workingDirectory, File sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
