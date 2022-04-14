package util.Configs.ConfigValidators;

import util.Configs.ConfigTypes.AbstractConfig;

public class IntegerRangeValidator implements Validator<Integer>
{
    private final int min, max;

    public IntegerRangeValidator(int min, int max)
    {
        this.min = min;
        this.max = max;
    }

    @Override public boolean validate(AbstractConfig<Integer> config)
    {
        if (config.get() == null)
        {
            return true;
        }
        return min <= config.get() && config.get() <= max;
    }

    @Override public String toString()
    {
        return "Value has to be between " + min + " and " + max + " (borders included)";
    }
}
