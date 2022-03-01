package util.Configs.Modules.Tepic;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class Tepic extends AbstractModule
{
    public FileStructure fileStructure;

    public final Config<File> inputDirectory = new Config<>(File.class);
    public final Config<File> inputPrevious = new Config<>(File.class);
    public final Config<File> inputOriginal = new Config<>(File.class);

    public final Config<File> inputReferenceGenome = new Config<>(File.class);
    public final Config<File> pathPwms = new Config<>(File.class);
    public final Config<Integer> cores = new Config<>(1);
    public final Config<File> bedChromatinSignal = new Config<>(File.class);
    public final Config<Integer> columnBedfile = new Config<>(Integer.class);
    public final Config<File> geneAnnotationFile = new Config<>(File.class);
    public final Config<Integer> windowSize = new Config<>(50000);
    public final Config<File> onlyDNasePeaks = new Config<>(File.class);
    public final Config<Boolean> exponentialDecay = new Config<>(true);
    public final Config<Boolean> doNotNormalizePeakLength = new Config<>(false);
    public final Config<Boolean> doNotGenerate = new Config<>(false);
    public final Config<Boolean> originalDecay = new Config<>(false);
    public final Config<File> psemsLengthFile = new Config<>(File.class);
    public final Config<Boolean> entireGeneBody = new Config<>(false);
    public final Config<Boolean> doZip = new Config<>(false);
    public final Config<Boolean> backgroundSequencesDirectory = new Config<>(false);
    public final Config<File> twoBitFile = new Config<>(File.class);
    public final Config<Double> pValue = new Config<>(0.05);
    public final Config<Integer> maxMinutesPerChromosome = new Config<>(3);
    public final Config<Boolean> chromosomePrefix = new Config<>(false);
    public final Config<Boolean> transcriptBased = new Config<>(false);
    public final Config<File> loopListFile = new Config<>(File.class);
    public final Config<Integer> loopWindows = new Config<>(5000);
    public final Config<Boolean> onlyPeakFeatures = new Config<>(false);
    public final Config<Integer> tpmCutoff = new Config<>(Integer.class);
    public final Config<File> ensgSymbolFile = new Config<>(File.class);
    public final Config<Boolean> tgeneTargetGenes = new Config<>(true);
    public final Config<Boolean> randomizeTfGeneMatrix = new Config<>(false);
    public final Config<String> tfBindingSiteSearch = new Config<>("INSIDE");
    public final Config<Integer> betweenMaxBps = new Config<>(500);

    public Tepic(File workingDirectory, File sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
