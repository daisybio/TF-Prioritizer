package util.Configs.Modules.FileStructure;

import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class FileStructure extends AbstractModule
{
    public MixOption mixOption;
    public BlacklistedRegions blacklistedRegions;
    public DeSeq2 deSeq2;
    public Tepic tepic;
    public Tgen tgen;
    public Dynamite dynamite;
    public DistributionAnalysis distributionAnalysis;
    public Igv igv;

    public FileStructure(File workingDirectory, File sourceDirectory, Logger logger)
            throws IllegalAccessException, InvocationTargetException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);

        init();
    }
}
