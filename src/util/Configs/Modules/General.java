package util.Configs.Modules;

import util.Configs.Config;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class General extends AbstractModule
{
    public final Config<Boolean> fileLogging = new Config<>(true, true);
    public final Config<File> logFile = extend(workingDirectory, "logfile.txt");
    public final Config<Boolean> ensgMappingEnabled = new Config<>(true, true);
    public final Config<Boolean> calculateTpmLengthsEnabled = new Config<>(true, true);
    public final Config<Boolean> calculateGenePositionsEnabled = new Config<>(true, true);

    public final Config<String> differentTps = new Config<>("DIFFERENT_TPS");
    public final Config<String> sameTps = new Config<>("SAME_TPS");

    public final Config<String> shebang = new Config<>("#!/bin/bash");

    public final Config<Boolean> performAllBackgroundDistributionAnalysis = new Config<>(false);
    public final Config<String> distributionAnalysisAllName = new Config<>("ALL");

    public final Config<Integer> threadLimit = new Config<>(10, true);

    public final Config<Boolean> redirectExternalScriptOutputStream = new Config<>(false, true);
    public final Config<Boolean> redirectExternalScriptErrorStream = new Config<>(false, true);

    public General(Config<File> workingDirectory, Config<File> sourceDirectory, Logger logger)
            throws IllegalAccessException, InvocationTargetException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
