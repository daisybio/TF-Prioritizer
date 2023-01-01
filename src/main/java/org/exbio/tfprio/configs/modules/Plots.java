package org.exbio.tfprio.configs.modules;

import org.exbio.pipejar.configs.ClassGetter;
import org.exbio.pipejar.configs.ConfigModule;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.ExternalConfig;

import java.util.List;

public class Plots extends ConfigModule {
    public final ExternalConfig<List<Double>> thresholds = new ExternalConfig<>(ClassGetter.getList());
}
