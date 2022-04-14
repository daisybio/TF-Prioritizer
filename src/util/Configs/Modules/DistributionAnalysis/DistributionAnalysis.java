package util.Configs.Modules.DistributionAnalysis;

import util.Configs.ConfigTypes.*;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class DistributionAnalysis extends AbstractModule
{
    public FileStructure fileStructure;

    public final InputConfig<Boolean> performAllAnalysis = new InputConfig<>(Boolean.class);
    public final InternalConfig<String> allName = new InternalConfig<>("ALL");

    public DistributionAnalysis(GeneratedFileStructure workingDirectory, SourceDirectoryFileStructure sourceDirectory,
                                Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
