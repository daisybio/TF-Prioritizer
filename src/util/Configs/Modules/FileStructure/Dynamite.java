package util.Configs.Modules.FileStructure;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class Dynamite extends AbstractModule
{
    public final Config<File> preprocessing = extend(workingDirectory, "05_A_DYNAMITE_preprocessing");
    public final Config<File> preprocessing_integrateData = extend(preprocessing, "integrateData");
    public final Config<File> preprocessing_prepareClassification = extend(preprocessing, "prepareClassification");
    public final Config<File> preprocessing_installRequiredPackages =
            extend(preprocessing, "X_install_required_packages");

    public final Config<File> output = extend(workingDirectory, "05_B_DYNAMITE_output");

    public Dynamite(File workingDirectory, File sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
