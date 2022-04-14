package util.Configs.ConfigTypes;

import util.FileManagement;
import util.Logger;

import java.io.File;
import java.io.IOException;

public class GeneratedFileStructure extends AbstractConfig<File>
{
    public GeneratedFileStructure(String path)
    {
        this(new File(path));
    }

    public GeneratedFileStructure(File file)
    {
        setValue(file);
    }

    private String noGenerationReason = null;

    @Override public boolean isWriteable()
    {
        return false;
    }

    @Override public boolean isValid(Logger logger)
    {
        return true;
    }

    @Override public boolean isSet()
    {
        return true;
    }

    public String getNoGenerationReason()
    {
        return noGenerationReason;
    }

    public void deleteAndSetNoGenerationReason(String noGenerationReason)
    {
        tryDelete();
        this.noGenerationReason = noGenerationReason;
    }

    @Override public void setValueObject(Object valueObject) throws IllegalAccessException
    {
        throw new IllegalAccessException("Trying to change a generated file structure via config file: " + name);
    }

    private void tryDelete()
    {
        try
        {
            FileManagement.deleteFileStructure(get());
        } catch (IOException ignore)
        {
        }
    }
}
