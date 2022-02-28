package util.Configs.Modules.FileStructure;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class MixOption extends AbstractModule
{
    public final Config<File> root = extend(workingDirectory, "00_MIX_OPTION");
    public final Config<File> preprocessingCheckChr = extend(root, "00_chromosome_annotation_checking");
    public final Config<File> sampleMixPreprocessing = extend(root, "00_PREPROCESSING_SAMPLE_MIX");
    public final Config<File> sampleMix = extend(root, "01_SAMPLE_MIX");
    public final Config<File> preprocessingHmMix = extend(root, "02_PREPROCESS_HM_MIX");
    public final Config<File> hmMix = extend(root, "03_HM_MIX");
    public final Config<File> footprintsBetweenPeaks = extend(root, "04_FOOTPRINTS");
    public final Config<File> mutuallyExclusive = extend(root, "02_MUTUALLY_EXCLUSIVE");
    public final Config<File> mutuallyExclusive_preprocessing = extend(mutuallyExclusive, "01_PREPROCESSING");
    public final Config<File> mutuallyExclusive_input = extend(mutuallyExclusive, "02_NEW_INPUT");

    public MixOption(File workingDirectory, File sourceDirectory, Logger logger)
            throws IllegalAccessException, InvocationTargetException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
