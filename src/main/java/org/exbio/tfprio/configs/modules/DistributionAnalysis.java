package org.exbio.tfprio.configs.modules;

import org.exbio.pipejar.configs.ConfigModule;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.ExternalConfig;

public class DistributionAnalysis extends ConfigModule {
    public final ExternalConfig<Boolean> scoreIncludeCounts = new ExternalConfig<>(Boolean.class);
}
