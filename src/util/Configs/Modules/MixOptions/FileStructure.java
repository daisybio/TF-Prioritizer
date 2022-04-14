package util.Configs.Modules.MixOptions;

import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.Configs.ConfigTypes.InputFileStructure;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class FileStructure extends AbstractModule
{
    public final GeneratedFileStructure d_root = extend(workingDirectory, "00_MIX_OPTION");
    public final GeneratedFileStructure d_preprocessingCheckChr = extend(d_root, "00_chromosome_annotation_checking");
    public final GeneratedFileStructure d_sampleMixPreprocessing = extend(d_root, "00_PREPROCESSING_SAMPLE_MIX");
    public final GeneratedFileStructure d_sampleMix = extend(d_root, "01_SAMPLE_MIX");
    public final GeneratedFileStructure d_preprocessingHmMix = extend(d_root, "02_PREPROCESS_HM_MIX");
    public final GeneratedFileStructure d_hmMix = extend(d_root, "03_HM_MIX");
    public final GeneratedFileStructure d_footprintsBetweenPeaks = extend(d_root, "04_FOOTPRINTS");
    public final GeneratedFileStructure d_mutuallyExclusive = extend(d_root, "02_MUTUALLY_EXCLUSIVE");
    public final GeneratedFileStructure d_mutuallyExclusive_preprocessing =
            extend(d_mutuallyExclusive, "01_PREPROCESSING");
    public final GeneratedFileStructure d_mutuallyExclusive_input = extend(d_mutuallyExclusive, "02_NEW_INPUT");

    public FileStructure(GeneratedFileStructure workingDirectory, InputFileStructure sourceDirectory, Logger logger)
            throws IllegalAccessException, InvocationTargetException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
