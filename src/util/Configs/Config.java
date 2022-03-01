package util.Configs;

import org.json.JSONObject;

import java.io.File;

public class Config<Type>
{
    private final Type defaultValue;
    private Type actualValue;
    private final boolean writeable;

    public Config()
    {
        this(null);
    }

    public Config(Type defaultValue)
    {
        this(defaultValue, false);
    }

    public Config(Type defaultValue, boolean writeable)
    {
        this.actualValue = this.defaultValue = defaultValue;
        this.writeable = writeable;
    }

    public void setValue(Object value) throws IllegalAccessException, ClassCastException
    {
        if (writeable)
        {
            if (value.getClass().equals(actualValue.getClass()))
            {
                try
                {
                    this.actualValue = (Type) value;
                } catch (Exception ignore)
                {
                }
            } else if (actualValue.getClass().equals(File.class) && value.getClass().equals(String.class))
            {
                actualValue = (Type) new File((String) value);
            } else
            {
                throw new ClassCastException(
                        "Trying to set a value of type " + value.getClass() + " to a config with" + " type " +
                                actualValue.getClass());
            }
        } else
        {
            throw new IllegalAccessException("Trying to manipulate an internal config.");
        }
    }

    public Type get()
    {
        return actualValue;
    }

    public String toString()
    {
        if (actualValue != null)
        {
            return actualValue.toString();
        }
        return "{NULL}";
    }

    public Object toJSONifyAble()
    {
        if (actualValue == null)
        {
            return JSONObject.NULL;
        }
        if (actualValue.getClass().equals(File.class))
        {
            return ((File) actualValue).getAbsolutePath();
        }
        return actualValue;
    }
}
