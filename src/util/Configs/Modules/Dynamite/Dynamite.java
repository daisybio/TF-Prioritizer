package util.Configs.Modules.Dynamite;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class Dynamite extends AbstractModule
{
    public FileStructure fileStructure;

    public final Config<Integer> preprocessing_IntegrateDataGeneIds = new Config<>(0);
    public final Config<Integer> preprocessing_IntegrateDataLog2fc = new Config<>(1);
    public final Config<File> preprocessing_IntegrateDataConsiderGeneFile = new Config<>(File.class);

    public final Config<String> outVar = new Config<>(String.class);
    public final Config<Integer> cores = new Config<>(1);
    public final Config<Double> alpha = new Config<>(0.1);
    public final Config<Double> testSize = new Config<>(0.2);
    public final Config<Integer> oFolds = new Config<>(3);
    public final Config<Integer> iFolds = new Config<>(6);
    public final Config<Boolean> balanced = new Config<>(true);
    public final Config<Boolean> performance = new Config<>(true);
    public final Config<Boolean> randomize = new Config<>(false);

    public Dynamite(Config<File> workingDirectory, Config<File> sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
