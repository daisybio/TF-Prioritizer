package util.Configs.Modules.Dynamite;

import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.Configs.ConfigTypes.InputFileStructure;
import util.Configs.ConfigTypes.InternalConfig;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class FileStructure extends AbstractModule
{
    public final GeneratedFileStructure d_preprocessing = extend(workingDirectory, "05_A_DYNAMITE_preprocessing");
    public final GeneratedFileStructure d_preprocessing_integrateData = extend(d_preprocessing, "integrateData");
    public final InternalConfig<String> s_preprocessing_integrateData_log2coeff =
            new InternalConfig<>("Integrated_Data_Log2_Quotient.txt");

    public final GeneratedFileStructure d_preprocessing_prepareClassification =
            extend(d_preprocessing, "prepareClassification");
    public final InternalConfig<String> s_preprocessing_prepareClassification_data =
            new InternalConfig<>("Integrated_Data_For_Classification.txt");

    public final GeneratedFileStructure d_preprocessing_installRequiredPackages =
            extend(d_preprocessing, "X_install_required_packages");
    public final GeneratedFileStructure s_preprocessing_installRequiredPackages_scrip =
            extend(d_preprocessing_installRequiredPackages, "install_required_packages_DYNAMITE.R");

    public final GeneratedFileStructure d_output = extend(workingDirectory, "05_B_DYNAMITE_output");
    public final InternalConfig<String> s_output_toBePlotted =
            new InternalConfig<>("Regression_Coefficients_Entire_Data_Set_Integrated_Data_For_Classification.txt");

    public FileStructure(GeneratedFileStructure workingDirectory, InputFileStructure sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
