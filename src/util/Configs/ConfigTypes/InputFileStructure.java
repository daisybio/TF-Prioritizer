package util.Configs.ConfigTypes;

import org.json.JSONObject;
import util.Logger;

import java.io.File;
import java.io.IOException;

public class InputFileStructure extends AbstractConfig<File>
{
    public InputFileStructure()
    {

    }

    public InputFileStructure(File file)
    {
        setValue(file);
    }

    @Override public boolean isWriteable()
    {
        return true;
    }

    @Override public boolean isValid(Logger logger)
    {
        return true;
    }

    @Override public boolean isSet()
    {
        return get() != null;
    }

    @Override public void setValueObject(Object valueObject) throws ClassCastException, IOException
    {
        if (valueObject != JSONObject.NULL)
        {
            if (valueObject.getClass().equals(String.class))
            {
                String valueString = (String) valueObject;
                setValue(new File(valueString));
            } else if (valueObject.getClass().equals(File.class))
            {
                setValue((File) valueObject);
            } else
            {
                throw new ClassCastException(
                        "Unknown data type: " + valueObject.getClass().getCanonicalName() + " for config: " + name);
            }

            if (!get().exists() || !get().canRead())
            {
                throw new IOException(
                        "Input file structure does not exist: " + get().getAbsolutePath() + " for config: " + name);
            }
        }
    }
}
