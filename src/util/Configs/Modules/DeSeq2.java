package util.Configs.Modules;

import util.Configs.Config;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class DeSeq2 extends AbstractModule
{
    public final Config<File> inputDirectory = new Config<>();
    public final Config<String> inputGeneID = new Config<>();
    public final Config<String> biomartDatasetSpecies = new Config<>();
    public final Config<String> biomartDatasetSymbolColumn = new Config<>();
    public final Config<Integer> countThreshold = new Config<>();
    public final Config<Double> tpmFilter = new Config<>();

    public DeSeq2(File workingDirectory, File sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
