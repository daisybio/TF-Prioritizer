package util.Configs.ConfigValidators;

import util.Configs.ConfigTypes.AbstractConfig;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class StringListValidator implements Validator<List<String>>
{
    private final Set<String> allowedValues;

    public StringListValidator(String... allowedValues)
    {
        this.allowedValues = new HashSet<>(List.of(allowedValues));
    }

    @Override public boolean validate(AbstractConfig<List<String>> config)
    {
        return allowedValues.containsAll(config.get());
    }

    @Override public String toString()
    {
        return "All values have to be member of this list: [" + String.join(", ", allowedValues) + "]";
    }
}
