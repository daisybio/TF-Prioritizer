package util.Configs.Modules.DeSeq2;

import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.Configs.ConfigTypes.InputConfig;
import util.Configs.ConfigTypes.InputFileStructure;
import util.Configs.ConfigTypes.SourceDirectoryFileStructure;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class DeSeq2 extends AbstractModule
{
    public FileStructure fileStructure;

    public final InputFileStructure inputDirectory = new InputFileStructure();
    public final InputFileStructure inputGeneID = new InputFileStructure();
    public final InputConfig<String> biomartDatasetSpecies = new InputConfig<>(String.class);
    public final InputConfig<String> biomartDatasetSymbolColumn = new InputConfig<>(String.class);
    public final InputConfig<Integer> countThreshold = new InputConfig<>(Integer.class);
    public final InputConfig<Double> tpmFilter = new InputConfig<>(Double.class);
    public final InputFileStructure d_enhancerDB = extend(extDirectory, "Enhancers_DB");

    public DeSeq2(GeneratedFileStructure workingDirectory, SourceDirectoryFileStructure sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
