package util.Configs.Modules.DeSeq2;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class DeSeq2 extends AbstractModule
{
    public FileStructure fileStructure;

    public final Config<File> inputDirectory = new Config<>(File.class);
    public final Config<String> inputGeneID = new Config<>(String.class);
    public final Config<String> biomartDatasetSpecies = new Config<>(String.class);
    public final Config<String> biomartDatasetSymbolColumn = new Config<>(String.class);
    public final Config<Integer> countThreshold = new Config<>(Integer.class);
    public final Config<Double> tpmFilter = new Config<>(Double.class);

    public DeSeq2(File workingDirectory, File sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
