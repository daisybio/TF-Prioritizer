package util.Configs.Modules.Igv;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static util.FileManagement.extend;

public class Igv extends AbstractModule
{
    public FileStructure fileStructure;

    public final Config<File> pathToIGV = extend(sourceDirectory, "IGV_2.11.2");
    public final Config<Integer> topLog2fc = new Config<>(15, true);
    public final Config<Boolean> topLog2fcIncludeLncRnaPseudogenes = new Config<>(true, true);
    public final Config<List> includePredictionData = new Config<>(List.class);
    public final Config<List> importantLociAllPrioTf = new Config<>(List.class);
    public final Config<File> pathToTfChipSeq = new Config<>(File.class);
    public final Config<File> pathToTdf = new Config<>(File.class);
    public final Config<List> enhancerDatabases = new Config<>(List.class);
    public final Config<Map> grcSynonymDict = new Config<>(Map.class);
    public final Config<Integer> port = new Config<>(60151, true);
    public final Config<String> speciesReferenceGenome = new Config<>(String.class);

    public Igv(Config<File> workingDirectory, Config<File> sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
