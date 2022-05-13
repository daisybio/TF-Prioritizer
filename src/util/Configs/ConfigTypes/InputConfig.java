package util.Configs.ConfigTypes;

import org.json.JSONArray;
import org.json.JSONObject;
import util.Configs.ConfigValidators.Validator;
import util.Logger;

import java.math.BigDecimal;
import java.util.*;

public class InputConfig<T> extends AbstractConfig<T>
{
    Set<Validator<T>> validators = new HashSet<>();

    private final Class<?> configClass;

    @SafeVarargs public InputConfig(Class<? extends T> configClass, Validator<T>... validators)
    {
        this.configClass = configClass;
        this.validators.addAll(List.of(validators));
    }

    @Override public boolean isWriteable()
    {
        return true;
    }

    @Override public boolean isValid(Logger logger)
    {
        boolean allValidatorsPassed = true;

        for (Validator<T> validator : validators)
        {
            boolean passed = validator.validate(this);
            if (!passed)
            {
                logger.warn(name + " does not match requirements: " + validator);
            }
            allValidatorsPassed = allValidatorsPassed && passed;
        }

        return allValidatorsPassed;
    }

    @Override public boolean isSet()
    {
        return get() != null;
    }

    @Override public void setValueObject(Object value) throws ClassCastException
    {
        if (value.getClass().equals(configClass))
        {
            setValue((T) value);
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
        } else if (configClass.equals(Double.class) && value.getClass().equals(Integer.class))
        {
            setValue((T) value);
        } else if (value == JSONObject.NULL)
        {
            setValue(null);
        } else
        {
            throw new ClassCastException(
                    "Trying to set a value of type " + value.getClass() + " to a config with type " + configClass +
                            ": " + getName());
        }
    }
}
