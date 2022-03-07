package util.Configs;

import org.json.JSONArray;
import org.json.JSONObject;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.math.BigDecimal;

public class Config<T>
{
    private final T defaultValue;
    private T actualValue;
    private final boolean writeable;
    private final Class configClass;

    public Config(Class<T> configClass)
    {
        this.actualValue = this.defaultValue = null;
        this.writeable = true;
        this.configClass = configClass;
    }

    public Config(T defaultValue)
    {
        this(defaultValue, false);
    }

    public Config(T defaultValue, boolean writeable)
    {
        assert defaultValue != null;
        this.actualValue = this.defaultValue = defaultValue;
        this.writeable = writeable;
        this.configClass = defaultValue.getClass();
    }

    public void setValue(Object value) throws IllegalAccessException, ClassCastException
    {
        if (writeable)
        {
            if (value.getClass().equals(configClass))
            {
                try
                {
                    this.actualValue = (T) value;
                } catch (Exception ignore)
                {
                }
            } else if (configClass.equals(File.class) && value.getClass().equals(String.class))
            {
                actualValue = (T) new File((String) value);
            } else if ((configClass.equals(List.class) || configClass.equals(java.util.ArrayList.class)) &&
                    value.getClass().equals(JSONArray.class))
            {
                List<?> bigDecimalList = ((JSONArray) value).toList();
                if (bigDecimalList.size() > 0 && bigDecimalList.get(0).getClass().equals(BigDecimal.class))
                {
                    List<Double> doubleList = new ArrayList<>();
                    for (Object bigDecimalValue : bigDecimalList)
                    {
                        doubleList.add(Double.valueOf(((BigDecimal) bigDecimalValue).doubleValue()));
                    }
                    actualValue = (T) doubleList;
                } else
                {
                    actualValue = (T) bigDecimalList;
                }
            } else if (configClass.equals(Map.class) && value.getClass().equals(JSONObject.class))
            {
                actualValue = (T) ((JSONObject) value).toMap();
            } else if (configClass.equals(Double.class) && value.getClass().equals(BigDecimal.class))
            {
                actualValue = (T) Double.valueOf(((BigDecimal) value).doubleValue());
            } else if (value == JSONObject.NULL)
            {
                actualValue = null;
            } else
            {
                throw new ClassCastException(
                        "Trying to set a value of type " + value.getClass() + " to a config with type " + configClass);
            }
        } else
        {
            throw new IllegalAccessException("Trying to manipulate an internal config.");
        }
    }

    public T get()
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

    public boolean isWriteable()
    {
        return writeable;
    }

    public boolean isSet()
    {
        return actualValue != null;
    }
}
