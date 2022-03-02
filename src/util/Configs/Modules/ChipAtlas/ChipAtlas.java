package util.Configs.Modules.ChipAtlas;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class ChipAtlas extends AbstractModule
{
    public FileStructure fileStructure;

    public final Config<Boolean> isEnabled = new Config<>(false);
    public final Config<String> urlToList =
            new Config<>("http://togodb.biosciencedbc.jp/togodb/release/chip_atlas_file_list.csv");
    public final Config<String> column_GeneVersion = new Config<>(".*genome.*assembly.*");
    public final Config<String> column_AntigenClass = new Config<>(".*antigen.*class.*");
    public final Config<String> column_Antigen = new Config<>(".*antigen.*");
    public final Config<String> column_CellTypeClass = new Config<>(".*cell.*type.*class.*");
    public final Config<String> column_url = new Config<>(".*file.*url.*");
    public final Config<String> genomeVersion = new Config<>(String.class);
    public final Config<String> tissueType = new Config<>(String.class);

    public ChipAtlas(Config<File> workingDirectory, Config<File> sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
