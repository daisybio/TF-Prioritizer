package org.exbio.tfprio.configs.modules;

import org.exbio.pipejar.configs.ConfigModule;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.ExternalConfig;
import org.exbio.pipejar.configs.ConfigValidators.PositiveDoubleValidator;

public class DeSeq2 extends ConfigModule {
    public final ExternalConfig<String> species = new ExternalConfig<>(String.class);
    public final ExternalConfig<String> speciesRefGenome = new ExternalConfig<>(String.class);
    public final ExternalConfig<Double> tpmFilter = new ExternalConfig<>(Double.class, new PositiveDoubleValidator());
}
