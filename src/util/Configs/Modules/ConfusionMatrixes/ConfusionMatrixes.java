package util.Configs.Modules.ConfusionMatrixes;

import tfprio.tfprio.TFPRIO;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.Configs.ConfigTypes.SourceDirectoryFileStructure;
import util.Configs.Configs;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.lang.reflect.InvocationTargetException;

public class ConfusionMatrixes extends AbstractModule
{
    public FileStructure fileStructure;

    /**
     * The default constructor.
     *
     * @param workingDirectory the {@link TFPRIO} working directory
     * @param sourceDirectory  the {@link TFPRIO} source directory
     * @param logger           the {@link Configs} logger
     */
    public ConfusionMatrixes(GeneratedFileStructure workingDirectory, SourceDirectoryFileStructure sourceDirectory,
                             Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
