package util.Configs.ConfigTypes;

import java.io.File;

public class SourceDirectoryFileStructure extends InputFileStructure
{
    public SourceDirectoryFileStructure(File file)
    {
        super(file);
    }

    @Override public boolean isWriteable()
    {
        return false;
    }
}
