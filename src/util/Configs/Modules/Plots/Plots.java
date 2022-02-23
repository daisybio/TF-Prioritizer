package util.Configs.Modules.Plots;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;
import java.util.Arrays;
import java.util.List;

public class Plots extends AbstractModule
{
    public FileStructure fileStructure;

    public final Config<List<Double>> thresholds = new Config<>(Arrays.asList(0.1, 0.2, 0.3, 0.4, 0.5, 0.6));
    public final Config<Integer> cutoffTps = new Config<>(1);
    public final Config<Integer> cutoffHms = new Config<>(1);
    public final Config<Integer> cutoffGcs = new Config<>(100);
    public final Config<Double> cutoffTpms = new Config<>(0.0);
    public final Config<Integer> topKGenes = new Config<>(30);
    public final Config<Double> mannWhitneyUPvalueCutoff = new Config<>(0.01);
    public final Config<String> distributionAnalysisScoreType = new Config<>("ECVL_GENE_COUNTS");
    public final Config<Double> trapPredictedSequenceLogosAffinityCutoff = new Config<>(0.05);

    public Plots(Config<File> workingDirectory, Config<File> sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}