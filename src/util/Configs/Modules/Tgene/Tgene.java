package util.Configs.Modules.Tgene;

import util.Configs.ConfigTypes.*;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class Tgene extends AbstractModule
{
    public FileStructure fileStructure;

    public final InputFileStructure pathToExecutable = new InputFileStructure();

    public final InputConfig<Boolean> noClosestLocus = new InputConfig<>(Boolean.class);
    public final InputConfig<Boolean> noClosestTss = new InputConfig<>(Boolean.class);
    public final InputConfig<Integer> maxLinkDistance = new InputConfig<>(Integer.class);
    public final InputConfig<Double> pValue = new InputConfig<>(Double.class);
    public final InputConfig<Boolean> selfRegulatory = new InputConfig<>(Boolean.class);
    public final InputConfig<Double> consensus = new InputConfig<>(Double.class);
    public final InputConfig<String> consensusCalc = new InputConfig<>(String.class);
    public final InputConfig<String> mtWriting = new InputConfig<>(String.class);

    public Tgene(GeneratedFileStructure workingDirectory, SourceDirectoryFileStructure sourceDirectory, Logger logger)
            throws IllegalAccessException, InvocationTargetException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
