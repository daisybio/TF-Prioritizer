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
    protected final T defaultValue;
    protected T actualValue;
    protected final boolean writeable;
    protected final Class configClass;
    protected final List<T> acceptedValues;
    private String name;

    public Config(Class<T> configClass)
    {
        this(configClass, null);
    }

    public Config(Class<T> configClass, List<T> acceptedValues)
    {
        this.actualValue = this.defaultValue = null;
        this.writeable = true;
        this.configClass = configClass;
        this.acceptedValues = acceptedValues;
    }

    public Config(T defaultValue)
    {
        this(defaultValue, false);
    }

    public Config(T defaultValue, boolean writeable)
    {
        this(defaultValue, null, writeable);
    }

    public Config(T defaultValue, List<T> acceptedValues, boolean writeable)
    {
        assert defaultValue != null;
        this.actualValue = this.defaultValue = defaultValue;
        this.writeable = writeable;
        this.configClass = defaultValue.getClass();
        this.acceptedValues = acceptedValues;
    }

    private void setActualValue(T value) throws IllegalArgumentException
    {
        if (acceptedValues != null)
        {
            if (acceptedValues.contains(value))
            {
                actualValue = value;
            } else
            {
                throw new IllegalArgumentException(
                        "The value " + value.toString() + " is not contained in the " + "accepted values: " +
                                acceptedValues);
            }
        } else
        {
            actualValue = value;
        }
    }

    public void setValue(Object value) throws IllegalAccessException, ClassCastException
    {
        if (isWriteable())
        {
            if (value.getClass().equals(configClass))
            {
                setActualValue((T) value);
            } else if (configClass.equals(File.class) && value.getClass().equals(String.class))
            {
                setActualValue((T) new File((String) value));
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
                    setActualValue((T) doubleList);
                } else
                {
                    setActualValue((T) bigDecimalList);
                }
            } else if (configClass.equals(Map.class) && value.getClass().equals(JSONObject.class))
            {
                setActualValue((T) ((JSONObject) value).toMap());
            } else if (configClass.equals(Double.class) && value.getClass().equals(BigDecimal.class))
            {
                setActualValue((T) Double.valueOf(((BigDecimal) value).doubleValue()));
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

    public void setName(String name)
    {
        this.name = name;
    }

    public String getName()
    {
        return name;
    }

    public boolean isWriteable()
    {
        return writeable;
    }

    public boolean isSet()
    {
        return actualValue != null;
    }

    public boolean isValid()
    {
        return true;
    }
}
