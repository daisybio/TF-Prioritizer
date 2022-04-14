package util.Configs.Modules.Tepic;

import util.Configs.ConfigTypes.*;
import util.Configs.ConfigValidators.IntegerRangeValidator;
import util.Configs.ConfigValidators.StringValidator;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class Tepic extends AbstractModule
{
    public FileStructure fileStructure;

    public final InputFileStructure d_ext_tepic = extend(extDirectory, "TEPIC", "TEPIC");
    public final InputFileStructure executable = extend(d_ext_tepic, "Code", "TEPIC.sh");
    public final InputFileStructure d_dynamiteScripts =
            extend(d_ext_tepic, "MachineLearningPipelines", "DYNAMITE", "Scripts");
    public final InputFileStructure f_dynamite_integrateDate = extend(d_dynamiteScripts, "integrateData.py");
    public final InputFileStructure f_dynamite_prepareForClassification =
            extend(d_dynamiteScripts, "prepareForClassification.R");
    public final InputFileStructure f_dynamite = extend(d_dynamiteScripts, "DYNAMITE.R");

    public final InputFileStructure inputDirectory = new InputFileStructure();

    public final InputFileStructure inputReferenceGenome = new InputFileStructure();
    public final InputFileStructure pathPwms = new InputFileStructure();
    public final InputFileStructure bedChromatinSignal = new InputFileStructure();
    public final InputConfig<Integer> columnBedfile = new InputConfig<>(Integer.class);
    public final InputFileStructure geneAnnotationFile = new InputFileStructure();
    public final InputConfig<Integer> windowSize = new InputConfig<>(Integer.class);
    public final InputFileStructure onlyDNasePeaks = new InputFileStructure();
    public final InputConfig<Boolean> exponentialDecay = new InputConfig<>(Boolean.class);
    public final InputConfig<Boolean> doNotNormalizePeakLength = new InputConfig<>(Boolean.class);
    public final InputConfig<Boolean> doNotGenerate = new InputConfig<>(Boolean.class);
    public final InputConfig<Boolean> originalDecay = new InputConfig<>(Boolean.class);
    public final InputFileStructure psemsLengthFile = new InputFileStructure();
    public final InputConfig<Boolean> entireGeneBody = new InputConfig<>(Boolean.class);
    public final InputConfig<Boolean> doZip = new InputConfig<>(Boolean.class);
    public final InputFileStructure backgroundSequencesDirectory = new InputFileStructure();
    public final InputFileStructure twoBitFile = new InputFileStructure();
    public final InputConfig<Double> pValue = new InputConfig<>(Double.class);
    public final InputConfig<Integer> maxMinutesPerChromosome = new InputConfig<>(Integer.class);
    public final InputConfig<Boolean> chromosomePrefix = new InputConfig<>(Boolean.class);
    public final InputConfig<Boolean> transcriptBased = new InputConfig<>(Boolean.class);
    public final InputFileStructure loopListFile = new InputFileStructure();
    public final InputConfig<Integer> loopWindows = new InputConfig<>(Integer.class);
    public final InputConfig<Boolean> onlyPeakFeatures = new InputConfig<>(Boolean.class);
    public final InputConfig<Integer> tpmCutoff =
            new InputConfig<>(Integer.class, new IntegerRangeValidator(1, Integer.MAX_VALUE));
    public final InputFileStructure ensgSymbolFile = new InputFileStructure();
    public final InputConfig<Boolean> tgeneTargetGenes = new InputConfig<>(Boolean.class);
    public final InputConfig<Boolean> randomizeTfGeneMatrix = new InputConfig<>(Boolean.class);
    public final InputConfig<String> tfBindingSiteSearch =
            new InputConfig<>(String.class, new StringValidator("INSIDE", "BETWEEN", "EXCL_BETWEEN"));

    // TODO: Add optional map to set bindingSiteSearch for each histone modification
    public final InputConfig<Integer> betweenMaxBps = new InputConfig<>(Integer.class);

    public Tepic(GeneratedFileStructure workingDirectory, SourceDirectoryFileStructure sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
