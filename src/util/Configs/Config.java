package util.Configs;

import org.json.JSONArray;
import org.json.JSONObject;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.math.BigDecimal;

/**
 * A single generic config containing an object which has an effect on the pipeline execution.
 *
 * @param <T> the data type of the config
 */
public class Config<T>
{
    /**
     * The config default value, is shown
     */
    protected final T defaultValue;

    /**
     * The config actual value
     */
    protected T actualValue;

    /**
     * Defines if the config can be overwritten by merged json files
     */
    protected final boolean writeable;

    /**
     * The class object of the default/actual value data elements.
     * <p>
     * Has to be stored like this, because this information cannot be obtained from the generic type.
     */
    protected final Class<?> configClass;

    /**
     * A list of possible config values.
     * <p>
     * Should be filled if only certain values are allowed for the config.
     */
    protected final List<T> acceptedValues;

    /**
     * The name of the config. Is used for indicating the user which config is invalid, if one is.
     */
    private String name;

    /**
     * Constructor for configs without a default value and no defined accepted values.
     * <p>
     * Writeable since not default value exists.
     *
     * @param configClass the class object of the generic data type
     */
    public Config(Class<T> configClass)
    {
        this(configClass, null);
    }

    /**
     * Constructor for configs without a default value and defined accepted values.
     * <p>
     * Writeable since no default value exists.
     *
     * @param configClass    the class object of the generic data type
     * @param acceptedValues the list of accepted values
     */
    public Config(Class<T> configClass, List<T> acceptedValues)
    {
        this.actualValue = this.defaultValue = null;
        this.writeable = true;
        this.configClass = configClass;
        this.acceptedValues = acceptedValues;
    }

    /**
     * Constructor for configs with a default value and no defined accepted values.
     * <p>
     * Not writeable by default.
     *
     * @param defaultValue the default value of the config
     */
    public Config(T defaultValue)
    {
        this(defaultValue, false);
    }

    /**
     * Constructor for configs with a default value and no defined accepted values, with writable option.
     *
     * @param defaultValue the default value of the config
     * @param writeable    defines if merged json files should be able to overwrite the config
     */
    public Config(T defaultValue, boolean writeable)
    {
        this(defaultValue, null, writeable);
    }

    /**
     * Constructor for configs with a default value, accepted values and writeable option.
     *
     * @param defaultValue   the default value of the config
     * @param acceptedValues the list of accepted values
     * @param writeable      defines if merged json files should be able to overwrite the config
     */
    public Config(T defaultValue, List<T> acceptedValues, boolean writeable)
    {
        assert defaultValue != null;
        this.actualValue = this.defaultValue = defaultValue;
        this.writeable = writeable;
        this.configClass = defaultValue.getClass();
        this.acceptedValues = acceptedValues;
    }

    /**
     * Set a certain value to the config.
     *
     * @param value the value to be set
     * @throws IllegalArgumentException if the value is not part of the {@link #acceptedValues}. Not thrown if the
     *                                  {@link #acceptedValues} are null.
     */
    private void setValue(T value) throws IllegalArgumentException
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

    /**
     * Try setting a value of object type to the config.
     * <p>
     * This is used during json merging since incoming data types are not clear at that point.
     *
     * @param value the value to be set
     * @throws IllegalAccessException if the config is not writeable
     * @throws ClassCastException     if the incoming data type does not match the generic data type
     */
    public void setValueObject(Object value) throws IllegalAccessException, ClassCastException
    {
        if (isWriteable())
        {
            if (value.getClass().equals(configClass))
            {
                setValue((T) value);
            } else if (configClass.equals(File.class) && value.getClass().equals(String.class))
            {
                setValue((T) new File((String) value));
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
                    setValue((T) doubleList);
                } else
                {
                    setValue((T) bigDecimalList);
                }
            } else if (configClass.equals(Map.class) && value.getClass().equals(JSONObject.class))
            {
                setValue((T) ((JSONObject) value).toMap());
            } else if (configClass.equals(Double.class) && value.getClass().equals(BigDecimal.class))
            {
                setValue((T) Double.valueOf(((BigDecimal) value).doubleValue()));
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

    /**
     * Get the value of the config.
     *
     * @return the config value
     */
    public T get()
    {
        return actualValue;
    }

    /**
     * Pretty-print the config value
     *
     * @return the config value string
     */
    @Override public String toString()
    {
        if (actualValue != null)
        {
            return actualValue.toString();
        }
        return "{NULL}";
    }

    /**
     * Get an object based on the config value which can then be converted to a JSONObject
     * <p>
     * Is used during config export.
     *
     * @return the JSONObject-convertible object
     */
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

    /**
     * Set the pretty name of this config.
     * <p>
     * Necessary because java does not allow accessing variable names at runtime
     *
     * @param name the pretty name to be set
     */
    public void setName(String name)
    {
        this.name = name;
    }

    /**
     * Get the pretty name of this config.
     *
     * @return the pretty name
     */
    public String getName()
    {
        return name;
    }

    /**
     * Check if the config is writable.
     *
     * @return the writable state of the config
     */
    public boolean isWriteable()
    {
        return writeable;
    }

    /**
     * Check if the config has got a value assigned.
     *
     * @return the assignment state
     */
    public boolean isSet()
    {
        return actualValue != null;
    }

    /**
     * Check if the config value does match the requirements.
     *
     * @return the validation state
     */
    public boolean isValid()
    {
        return true;
    }
}
