package util.Configs.Modules;

import util.Configs.Config;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;
import java.util.List;

public class Misc extends AbstractModule
{
    public final Config<String> analysisTypes_distribution = new Config<>("DISTRIBUTION_ANALYSIS");
    public final Config<String> analysisTypes_regression = new Config<>("REGRESSION_COEFFICIENT");

    public final Config<String> conversionHelpers_conversionInput = new Config<>(String.class);
    public final Config<String> conversionHelpers_inputFormat = new Config<>(String.class);
    public final Config<List> conversionHelpers_groups = new Config<>(List.class);
    public final Config<Boolean> conversionHelpers_cutToInteger = new Config<>(false);
    public final Config<String> conversionHelpers_biomartColumn = new Config<>(String.class);

    public Misc(Config<File> workingDirectory, Config<File> sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
