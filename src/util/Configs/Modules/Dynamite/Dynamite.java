package util.Configs.Modules.Dynamite;

import util.Configs.ConfigTypes.*;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class Dynamite extends AbstractModule
{
    public FileStructure fileStructure;

    public final InternalConfig<Integer> preprocessing_IntegrateDataGeneIds = new InternalConfig<>(0);
    public final InternalConfig<Integer> preprocessing_IntegrateDataLog2fc = new InternalConfig<>(1);
    public final InputFileStructure preprocessing_IntegrateDataConsiderGeneFile = new InputFileStructure();

    public final InputConfig<String> outVar = new InputConfig<>(String.class);
    public final InputConfig<Integer> cores = new InputConfig<>(Integer.class);
    public final InputConfig<Double> alpha = new InputConfig<>(Double.class);
    public final InputConfig<Double> testSize = new InputConfig<>(Double.class);
    public final InputConfig<Integer> oFolds = new InputConfig<>(Integer.class);
    public final InputConfig<Integer> iFolds = new InputConfig<>(Integer.class);
    public final InputConfig<Boolean> balanced = new InputConfig<>(Boolean.class);
    public final InputConfig<Boolean> performance = new InputConfig<>(Boolean.class);
    public final InputConfig<Boolean> randomize = new InputConfig<>(Boolean.class);

    public Dynamite(GeneratedFileStructure workingDirectory, SourceDirectoryFileStructure sourceDirectory,
                    Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
