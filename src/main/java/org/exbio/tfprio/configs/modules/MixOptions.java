package org.exbio.tfprio.configs.modules;

import org.exbio.pipejar.configs.ConfigModule;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.ExternalConfig;
import org.exbio.pipejar.configs.ConfigValidators.RequiresOtherToBeTrueValidator;
import org.exbio.pipejar.configs.ConfigValidators.StringValidator;

import java.io.File;

public class MixOptions extends ConfigModule {
    public final ExternalConfig<Boolean> perform = new ExternalConfig<>(Boolean.class);
    public final ExternalConfig<Integer> minOccurrence = new ExternalConfig<>(Integer.class);
    public final ExternalConfig<File> blackListPath = new ExternalConfig<>(File.class);
    public final ExternalConfig<String> biomartDatasetSpecies = new ExternalConfig<>(String.class);
    public final ExternalConfig<String> tfBindingSiteSearch =
            new ExternalConfig<>(String.class, new StringValidator("BETWEEN", "EXCL_BETWEEN", "INSIDE"));
    public final ExternalConfig<String> seqType =
            new ExternalConfig<>(String.class, new StringValidator("chip-seq", "atac-seq", "dnase-seq"));
    public final ExternalConfig<Integer> maxSpaceBetweenPeaks = new ExternalConfig<>(Integer.class);
    public final ExternalConfig<Boolean> mixMutuallyExclusive =
            new ExternalConfig<>(Boolean.class, new RequiresOtherToBeTrueValidator(perform));
    public final ExternalConfig<Boolean> differentialPeakSignals =
            new ExternalConfig<>(Boolean.class, new RequiresOtherToBeTrueValidator(mixMutuallyExclusive));
}
