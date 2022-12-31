package org.exbio.tfprio.configs.modules;

import org.exbio.pipejar.configs.ConfigModule;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.ExternalConfig;
import org.exbio.pipejar.configs.ConfigValidators.PositiveDoubleValidator;
import org.exbio.pipejar.configs.ConfigValidators.PositiveIntegerValidator;

public class DYNAMITE extends ConfigModule {
    public final ExternalConfig<String> outVar = new ExternalConfig<>(String.class);
    public final ExternalConfig<Integer> oFolds = new ExternalConfig<>(Integer.class, new PositiveIntegerValidator());
    public final ExternalConfig<Integer> iFolds = new ExternalConfig<>(Integer.class, new PositiveIntegerValidator());
    public final ExternalConfig<Boolean> performance = new ExternalConfig<>(Boolean.class);
    public final ExternalConfig<Double> alpha = new ExternalConfig<>(Double.class, new PositiveDoubleValidator());
    public final ExternalConfig<Boolean> randomize = new ExternalConfig<>(Boolean.class);
}
