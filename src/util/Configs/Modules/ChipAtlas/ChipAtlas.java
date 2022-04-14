package util.Configs.Modules.ChipAtlas;

import util.Configs.ConfigTypes.*;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class ChipAtlas extends AbstractModule
{
    public FileStructure fileStructure;

    public final InputConfig<Boolean> isEnabled = new InputConfig<>(Boolean.class);
    public final InternalConfig<String> urlToList =
            new InternalConfig<>("http://togodb.biosciencedbc.jp/togodb/release/chip_atlas_file_list.csv");
    public final InternalConfig<String> column_GeneVersion = new InternalConfig<>(".*genome.*assembly.*");
    public final InternalConfig<String> column_AntigenClass = new InternalConfig<>(".*antigen.*class.*");
    public final InternalConfig<String> column_Antigen = new InternalConfig<>(".*antigen.*");
    public final InternalConfig<String> column_CellTypeClass = new InternalConfig<>(".*cell.*type.*class.*");
    public final InternalConfig<String> column_url = new InternalConfig<>(".*file.*url.*");
    public final InputConfig<String> genomeVersion = new InputConfig<>(String.class);
    public final InputConfig<String> tissueType = new InputConfig<>(String.class);

    public ChipAtlas(GeneratedFileStructure workingDirectory, SourceDirectoryFileStructure sourceDirectory,
                     Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
