package util.Configs.Modules.Tgene;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class Tgene extends AbstractModule
{
    public FileStructure fileStructure;

    public final Config<File> pathToExecutable = new Config<>(File.class);

    public final Config<Boolean> noClosestLocus = new Config<>(false);
    public final Config<Boolean> noClosestTss = new Config<>(false);
    public final Config<Integer> maxLinkDistance = new Config<>(500000);
    public final Config<Double> pValue = new Config<>(0.05);
    public final Config<Boolean> selfRegulatory = new Config<>(false);
    public final Config<Double> consensus = new Config<>(.5);
    public final Config<String> consensusCalc = new Config<>("INCREASE_TGENE_TFS");
    public final Config<String> mtWriting = new Config<>("MT");

    public Tgene(File workingDirectory, File sourceDirectory, Logger logger)
            throws IllegalAccessException, InvocationTargetException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
