package util.Configs.ConfigValidators;

import util.Configs.ConfigTypes.AbstractConfig;

public interface Validator<T>
{
    boolean validate(AbstractConfig<T> config);

    String toString();
}
