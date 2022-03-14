package util.Configs.Modules.DeSeq2;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Configs.RequiredConfig;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class DeSeq2 extends AbstractModule
{
    public FileStructure fileStructure;

    public final Config<File> inputDirectory = new Config<>(File.class);
    public final Config<String> inputGeneID = new Config<>(String.class);
    public final Config<String> biomartDatasetSpecies = new Config<>(String.class);
    public final Config<String> biomartDatasetSymbolColumn = new RequiredConfig<>(String.class);
    public final Config<Integer> countThreshold = new Config<>(Integer.class);
    public final Config<Double> tpmFilter = new Config<>(Double.class);
    public final Config<File> d_enhancerDB = extend(extDirectory, "Enhancers_DB");

    public DeSeq2(Config<File> workingDirectory, Config<File> sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
