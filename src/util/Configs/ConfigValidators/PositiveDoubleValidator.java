package util.Configs.ConfigValidators;

import util.Configs.ConfigTypes.AbstractConfig;

public class PositiveDoubleValidator implements Validator<Double>
{
    @Override public boolean validate(AbstractConfig<Double> config)
    {
        if (config.isSet())
        {
            return config.get() > 0;
        }
        return true;
    }
}
