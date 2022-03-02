package util.Configs.Modules.ScriptTemplates;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class ScriptTemplates extends AbstractModule
{
    public final Config<File> d_root = extend(sourceDirectory, "scripts");
    public final Config<File> f_heatmaps = extend(sourceDirectory, "heatmap.R");

    public ScriptTemplates(Config<File> workingDirectory, Config<File> sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
