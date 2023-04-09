package org.exbio.tfprio.configs.modules;

import org.exbio.pipejar.configs.ConfigModule;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.ExternalConfig;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.InternalConfig;

import java.io.File;


public class EHMM extends ConfigModule {
    public final ExternalConfig<File> bamDirectory = new ExternalConfig<>(File.class);
    public final InternalConfig<Integer> nStates = new InternalConfig<>(10);
    public final InternalConfig<Double> pseudoCount = new InternalConfig<>(1e-6);
    public final InternalConfig<Integer> nBins = new InternalConfig<>(100);
    public final InternalConfig<Integer> nThreads = new InternalConfig<>(25);
}
