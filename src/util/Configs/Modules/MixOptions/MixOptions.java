package util.Configs.Modules.MixOptions;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Arrays;

public class MixOptions extends AbstractModule
{
    public FileStructure fileStructure;

    public final Config<String> level =
            new Config<>(String.class, new ArrayList<>(Arrays.asList("SAMPLE_LEVEL", "HM_LEVEL")));
    public final Config<String> option =
            new Config<>(String.class, new ArrayList<>(Arrays.asList("UNION", "INTERSECTION")));
    public final Config<Integer> occurrenceIntersection = new Config<>(2, true);
    public final Config<Boolean> mutuallyExclusive = new Config<>(false, true);
    public final Config<Boolean> mutuallyExclusiveDifferentialPeakSignals = new Config<>(true, true);

    public MixOptions(Config<File> workingDirectory, Config<File> sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }

    @Override public boolean validate()
    {
        if (!super.validate())
        {
            return false;
        }

        if (level.isSet())
        {
            if (!option.isSet())
            {
                logger.error("If MixOption>level is set, also an option has to be set.");
                return false;
            }
        }

        if (mutuallyExclusive.get())
        {
            if (!level.isSet())
            {
                logger.error("If MixOption>mutuallyExclusive is set, also MixOption>level has to be set.");
                return false;
            }
        }

        return true;
    }
}
