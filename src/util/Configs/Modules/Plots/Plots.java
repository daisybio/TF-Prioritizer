package util.Configs.Modules.Plots;

import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.Configs.ConfigTypes.InputConfig;
import util.Configs.ConfigTypes.InputFileStructure;
import util.Configs.ClassGetter;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.lang.reflect.InvocationTargetException;
import java.util.List;

public class Plots extends AbstractModule
{
    public FileStructure fileStructure;

    public final InputConfig<List<Double>> thresholds = new InputConfig<>(ClassGetter.getDoubleList());
    public final InputConfig<Integer> cutoffTps = new InputConfig<>(Integer.class);
    public final InputConfig<Integer> cutoffHms = new InputConfig<>(Integer.class);
    public final InputConfig<Integer> cutoffGcs = new InputConfig<>(Integer.class);
    public final InputConfig<Double> cutoffTpms = new InputConfig<>(Double.class);
    public final InputConfig<Integer> topKGenes = new InputConfig<>(Integer.class);
    public final InputConfig<Double> mannWhitneyUPvalueCutoff = new InputConfig<>(Double.class);
    public final InputConfig<String> distributionAnalysisScoreType = new InputConfig<>(String.class);
    public final InputConfig<Double> trapPredictedSequenceLogosAffinityCutoff = new InputConfig<>(Double.class);

    public Plots(GeneratedFileStructure workingDirectory, InputFileStructure sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
