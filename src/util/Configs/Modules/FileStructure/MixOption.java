package util.Configs.Modules.FileStructure;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class MixOption extends AbstractModule
{
    public final Config<File> d_root = extend(workingDirectory, "00_MIX_OPTION");
    public final Config<File> d_preprocessingCheckChr = extend(d_root, "00_chromosome_annotation_checking");
    public final Config<File> d_sampleMixPreprocessing = extend(d_root, "00_PREPROCESSING_SAMPLE_MIX");
    public final Config<File> d_sampleMix = extend(d_root, "01_SAMPLE_MIX");
    public final Config<File> d_preprocessingHmMix = extend(d_root, "02_PREPROCESS_HM_MIX");
    public final Config<File> d_hmMix = extend(d_root, "03_HM_MIX");
    public final Config<File> d_footprintsBetweenPeaks = extend(d_root, "04_FOOTPRINTS");
    public final Config<File> d_mutuallyExclusive = extend(d_root, "02_MUTUALLY_EXCLUSIVE");
    public final Config<File> d_mutuallyExclusive_preprocessing = extend(d_mutuallyExclusive, "01_PREPROCESSING");
    public final Config<File> d_mutuallyExclusive_input = extend(d_mutuallyExclusive, "02_NEW_INPUT");

    public MixOption(File workingDirectory, File sourceDirectory, Logger logger)
            throws IllegalAccessException, InvocationTargetException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
