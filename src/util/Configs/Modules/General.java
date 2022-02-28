package util.Configs.Modules;

import util.Configs.Config;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class General extends AbstractModule
{
    public final Config<Boolean> fileLogging = new Config<>(true, true);
    public final Config<Boolean> ensgMappingEnabled = new Config<>(true, true);
    public final Config<Boolean> calculateTpmLengthsEnabled = new Config<>(true, true);
    public final Config<Boolean> calculateGenePositionsEnabled = new Config<>(true, true);

    public final Config<String> differentTps = new Config<>("DIFFERENT_TPS", false);
    public final Config<String> sameTps = new Config<>("SAME_TPS", false);

    public final Config<File> d_ext =
            new Config<>(new File(sourceDirectory.getAbsolutePath() + File.separator + "ext"), false);

    public final Config<String> shebang = new Config<>("#!/bin/bash", false);

    public final Config<Boolean> performAllBackgroundDistributionAnalysis = new Config<>(false, false);
    public final Config<String> distributionAnalysisAllName = new Config<>("ALL", false);


    public General(File workingDirectory, File sourceDirectory, Logger logger)
            throws IllegalAccessException, InvocationTargetException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
