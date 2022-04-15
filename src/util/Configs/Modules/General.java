package util.Configs.Modules;

import util.Configs.ConfigTypes.*;
import util.Configs.ConfigValidators.PositiveIntegerValidator;
import util.Logger;

import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class General extends AbstractModule
{
    public final InputConfig<Boolean> fileLogging = new InputConfig<>(Boolean.class);
    public final GeneratedFileStructure logFile = extend(workingDirectory, "logfile.txt");
    public final InputConfig<Boolean> calculateTpmLengthsEnabled = new InputConfig<>(Boolean.class);
    public final InputConfig<Boolean> calculateGenePositionsEnabled = new InputConfig<>(Boolean.class);

    public final InternalConfig<String> differentTps = new InternalConfig<>("DIFFERENT_TPS");
    public final InternalConfig<String> sameTps = new InternalConfig<>("SAME_TPS");

    public final InputConfig<Integer> threadLimit = new InputConfig<>(Integer.class, new PositiveIntegerValidator());

    public final InputConfig<Boolean> redirectExternalScriptOutputStream = new InputConfig<>(Boolean.class);
    public final InputConfig<Boolean> redirectExternalScriptErrorStream = new InputConfig<>(Boolean.class);

    public final GeneratedFileStructure d_workflowHashes = extend(workingDirectory, ".hashes");

    public General(GeneratedFileStructure workingDirectory, SourceDirectoryFileStructure sourceDirectory, Logger logger)
            throws IllegalAccessException, InvocationTargetException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
