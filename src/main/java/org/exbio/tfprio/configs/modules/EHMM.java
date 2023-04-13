package org.exbio.tfprio.configs.modules;

import org.exbio.pipejar.configs.ConfigModule;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.ExternalConfig;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.InternalConfig;
import org.exbio.pipejar.configs.ConfigValidators.StringValidator;

import java.io.File;


public class EHMM extends ConfigModule {
    public final ExternalConfig<File> bamDirectory = new ExternalConfig<>(File.class);
    public final InternalConfig<String> experiment = new InternalConfig<>("chip-seq");
    public final InternalConfig<String> chromatinAccessClass = new InternalConfig<>("RNA polymerase");
    public final InternalConfig<String> chromatinAccessAntigen = new InternalConfig<>("RNA polymerase II");
    public final InternalConfig<Integer> nStates = new InternalConfig<>(10);
    public final InternalConfig<Double> pseudoCount = new InternalConfig<>(1.0);
    public final InternalConfig<Integer> nBins = new InternalConfig<>(10);
    public final InternalConfig<Integer> nThreads = new InternalConfig<>(25);
}
