package org.exbio.tfprio.configs.modules;

import org.exbio.pipejar.configs.ClassGetter;
import org.exbio.pipejar.configs.ConfigModule;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.ExternalConfig;
import org.exbio.pipejar.configs.ConfigValidators.PositiveIntegerValidator;

import java.util.List;

public class Plots extends ConfigModule {
    public final ExternalConfig<List<Double>> thresholds = new ExternalConfig<>(ClassGetter.getList());
    public final ExternalConfig<Integer> cutoffGroups = new ExternalConfig<>(Integer.class);
    public final ExternalConfig<Double> cutoffGroupCounts = new ExternalConfig<>(Double.class);
    public final ExternalConfig<Double> cutoffTPM = new ExternalConfig<>(Double.class);
    public final ExternalConfig<Integer> cutoffHms = new ExternalConfig<>(Integer.class);
    public final ExternalConfig<Integer> topKTargetGenes =
            new ExternalConfig<>(Integer.class, new PositiveIntegerValidator());
}
