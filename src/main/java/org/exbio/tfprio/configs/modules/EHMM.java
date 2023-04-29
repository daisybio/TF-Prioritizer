package org.exbio.tfprio.configs.modules;

import org.exbio.pipejar.configs.ConfigModule;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.ExternalConfig;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.InputConfig;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.InternalConfig;
import org.exbio.pipejar.configs.ConfigValidators.IntegerRangeValidator;
import org.exbio.pipejar.configs.ConfigValidators.SingleValidator;

import java.io.File;
import java.util.Set;

class DoubleValidator extends SingleValidator<Double>{
    Double min;
    Double max;

    public DoubleValidator(Double min, Double max) {
        this.min = min;
        this.max = max;
    }

    @Override
    public boolean validateSingle(InputConfig<Double> inputConfig) {
        if (inputConfig.get() == null) return true;
        return min < inputConfig.get() && inputConfig.get() < max;
    }
    @Override
    public String toString() {
        return "Value has to be between " + this.min + " and " + this.max + " (borders excluded)";
    }
}


public class EHMM extends ConfigModule {
    public final ExternalConfig<File> bamDirectory = new ExternalConfig<>(File.class);
    public final ExternalConfig<Integer> nStates = new ExternalConfig<>(Integer.class, new IntegerRangeValidator(3, 12));
    public final ExternalConfig<Integer> nSamples = new ExternalConfig<>(Integer.class, new IntegerRangeValidator(100, 100000));
    public final ExternalConfig<Double> trainSplit = new ExternalConfig<>(Double.class, new DoubleValidator(0d, 1d));
    public final ExternalConfig<Integer> topQuantile =
            new ExternalConfig<>(Integer.class, new IntegerRangeValidator(1, 4));
    public final ExternalConfig<Double> score = new ExternalConfig<>(Double.class, new DoubleValidator(0d, 1d));
    // Internal configs
    public final InternalConfig<Double> pseudoCount = new InternalConfig<>(1.0);
    public final InternalConfig<Set<String>> antigenClassKeys = new InternalConfig<>(
            Set.of("ATAC-Seq", "Histone_H3K27ac", "Histone_H3K4me1", "Histone_H3K4me3"));
    public final InternalConfig<String> threshold = new InternalConfig<>("5");
    public final InternalConfig<Integer> nBins = new InternalConfig<>(100);
    public final InternalConfig<Integer> nThreads = new InternalConfig<>(21);
    public final InternalConfig<Integer> genomicRegionSize = new InternalConfig<>(2000);
}
