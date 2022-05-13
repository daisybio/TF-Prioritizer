package util.Configs.Modules.AngularReport;

import tfprio.tfprio.TFPRIO;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.Configs.ConfigTypes.SourceDirectoryFileStructure;
import util.Configs.Configs;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class FileStructure extends AbstractModule
{
    public final GeneratedFileStructure d_preprocessing = extend(workingDirectory, "Report_Preprocessing");
    public final GeneratedFileStructure f_data = extend(d_preprocessing, "data.json");
    public final GeneratedFileStructure d_data = extend(d_preprocessing, "data");
    public final GeneratedFileStructure f_script = extend(d_preprocessing, "generate.sh");

    public final GeneratedFileStructure d_output = extend(workingDirectory, "AngularReport");

    /**
     * The default constructor.
     *
     * @param workingDirectory the {@link TFPRIO} working directory
     * @param sourceDirectory  the {@link TFPRIO} source directory
     * @param logger           the {@link Configs} logger
     */
    public FileStructure(GeneratedFileStructure workingDirectory, SourceDirectoryFileStructure sourceDirectory,
                         Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
