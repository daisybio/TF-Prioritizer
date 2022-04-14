package util.Configs.ConfigValidators;

import util.Configs.ConfigTypes.AbstractConfig;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class StringValidator implements Validator<String>
{
    private final Set<String> allowedValues;

    public StringValidator(String... allowedValues)
    {
        this.allowedValues = new HashSet<>(List.of(allowedValues));
    }

    @Override public boolean validate(AbstractConfig<String> config)
    {
        if (!config.isSet())
        {
            return true;
        }

        Object value = config.get();
        return allowedValues.contains(value);
    }

    @Override public String toString()
    {
        return "The following values are allowed: [" + String.join(", ", allowedValues) + "]";
    }
}
