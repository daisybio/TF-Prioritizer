package util.Configs.Modules.Blacklist;

import util.Configs.ClassGetter;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.Configs.ConfigTypes.InputConfig;
import util.Configs.ConfigTypes.InputFileStructure;
import util.Configs.ConfigValidators.StringListValidator;
import util.Configs.ConfigValidators.StringValidator;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.lang.reflect.InvocationTargetException;
import java.util.List;

public class Blacklist extends AbstractModule
{
    public FileStructure fileStructure;

    public final InputFileStructure bedFilePath = new InputFileStructure();
    public final InputConfig<List<String>> signalsToIgnore = new InputConfig<>(ClassGetter.getStringList(),
            new StringListValidator("Low_Mappability".toUpperCase(), "High_Signal_Region".toUpperCase()));

    public Blacklist(GeneratedFileStructure workingDirectory, InputFileStructure sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
