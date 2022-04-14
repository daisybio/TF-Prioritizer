package util.Configs.Modules.MixOptions;

import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.Configs.ConfigTypes.InputConfig;
import util.Configs.ConfigTypes.InputFileStructure;
import util.Configs.ConfigValidators.StringValidator;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Arrays;

public class MixOptions extends AbstractModule
{
    public FileStructure fileStructure;

    public final InputConfig<String> level =
            new InputConfig<>(String.class, new StringValidator("SAMPLE_LEVEL", "HM_LEVEL"));
    public final InputConfig<String> option =
            new InputConfig<>(String.class, new StringValidator("UNION", "INTERSECTION"));
    public final InputConfig<Integer> occurrenceIntersection = new InputConfig<>(Integer.class);
    public final InputConfig<Boolean> mutuallyExclusive = new InputConfig<>(Boolean.class);
    public final InputConfig<Boolean> mutuallyExclusiveDifferentialPeakSignals = new InputConfig<>(Boolean.class);

    public MixOptions(GeneratedFileStructure workingDirectory, InputFileStructure sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
