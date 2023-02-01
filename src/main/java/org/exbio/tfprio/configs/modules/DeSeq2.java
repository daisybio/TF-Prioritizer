package org.exbio.tfprio.configs.modules;

import org.exbio.pipejar.configs.ClassGetter;
import org.exbio.pipejar.configs.ConfigModule;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.ExternalConfig;
import org.exbio.pipejar.configs.ConfigValidators.PositiveDoubleValidator;
import org.exbio.pipejar.configs.ConfigValidators.PositiveIntegerValidator;

import java.util.Map;

public class DeSeq2 extends ConfigModule {
    public final ExternalConfig<Double> tpmFilter = new ExternalConfig<>(Double.class, new PositiveDoubleValidator());
    public final ExternalConfig<Integer> countFilter =
            new ExternalConfig<>(Integer.class, new PositiveIntegerValidator());
    public final ExternalConfig<Map<String, Integer>> batches = new ExternalConfig<>(ClassGetter.getMap());
}
