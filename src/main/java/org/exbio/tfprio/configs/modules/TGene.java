package org.exbio.tfprio.configs.modules;

import org.exbio.pipejar.configs.ConfigModule;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.ExternalConfig;

import java.io.File;

public class TGene extends ConfigModule {
    public final ExternalConfig<File> executable = new ExternalConfig<>(File.class);
    public final ExternalConfig<Boolean> noClosestLocus = new ExternalConfig<>(Boolean.class);
    public final ExternalConfig<Boolean> noClosestTss = new ExternalConfig<>(Boolean.class);
    public final ExternalConfig<Integer> maxLinkDistance = new ExternalConfig<>(Integer.class);
    public final ExternalConfig<Double> pValue = new ExternalConfig<>(Double.class);
}

