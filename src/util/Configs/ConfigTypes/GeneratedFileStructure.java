package util.Configs.ConfigTypes;

import util.Logger;

import java.io.File;

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

    private String noGenerationReason = "[no known reason]";

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

    public void setNoGenerationReason(String noGenerationReason)
    {
        this.noGenerationReason = noGenerationReason;
    }

    @Override public void setValueObject(Object valueObject) throws IllegalAccessException
    {
        throw new IllegalAccessException("Trying to change a generated file structure via config file: " + name);
    }
}
