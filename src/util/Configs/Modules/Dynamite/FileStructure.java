package util.Configs.Modules.Dynamite;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class FileStructure extends AbstractModule
{
    public final Config<File> d_preprocessing = extend(workingDirectory, "05_A_DYNAMITE_preprocessing");
    public final Config<File> d_preprocessing_integrateData = extend(d_preprocessing, "integrateData");
    public final Config<String> s_preprocessing_integrateData_log2coeff =
            new Config<>("Integrated_Data_Log2_Quotient.txt");

    public final Config<File> d_preprocessing_prepareClassification = extend(d_preprocessing, "prepareClassification");
    public final Config<String> s_preprocessing_prepareClassification_data =
            new Config<>("Integrated_Data_For_Classification.txt");

    public final Config<File> d_preprocessing_installRequiredPackages =
            extend(d_preprocessing, "X_install_required_packages");
    public final Config<File> s_preprocessing_installRequiredPackages_scrip =
            extend(d_preprocessing_installRequiredPackages, "install_required_packages_DYNAMITE.R");

    public final Config<File> d_output = extend(workingDirectory, "05_B_DYNAMITE_output");
    public final Config<String> s_output_toBePlotted =
            new Config<>("Regression_Coefficients_Entire_Data_Set_Integrated_Data_For_Classification.txt");

    public FileStructure(Config<File> workingDirectory, Config<File> sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
