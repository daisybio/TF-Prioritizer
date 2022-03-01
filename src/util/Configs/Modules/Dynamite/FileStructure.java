package util.Configs.Modules.Dynamite;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class FileStructure extends AbstractModule
{
    public final Config<File> d_preprocessing = extend(workingDirectory, "05_A_DYNAMITE_preprocessing");
    public final Config<File> d_preprocessing_integrateData = extend(d_preprocessing, "integrateData");
    public final Config<File> d_preprocessing_prepareClassification = extend(d_preprocessing, "prepareClassification");
    public final Config<File> d_preprocessing_installRequiredPackages =
            extend(d_preprocessing, "X_install_required_packages");

    public final Config<File> d_output = extend(workingDirectory, "05_B_DYNAMITE_output");

    public FileStructure(File workingDirectory, File sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
