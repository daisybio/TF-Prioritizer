package util.Configs.Modules.AngularReport;

import tfprio.tfprio.TFPRIO;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.Configs.ConfigTypes.SourceDirectoryFileStructure;
import util.Configs.Configs;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class AngularReport extends AbstractModule
{
    public FileStructure fileStructure;

    public final SourceDirectoryFileStructure d_angularSource = extend(sourceDirectory, "ext", "Report");

    /**
     * The default constructor.
     *
     * @param workingDirectory the {@link TFPRIO#workingDirectory}
     * @param sourceDirectory  the {@link TFPRIO#sourceDirectory}
     * @param logger           the {@link Configs} logger
     */
    public AngularReport(GeneratedFileStructure workingDirectory, SourceDirectoryFileStructure sourceDirectory,
                         Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
