package org.exbio.tfprio.configs.modules;

import org.exbio.pipejar.configs.ConfigModule;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.ExternalConfig;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.InternalConfig;

import java.io.File;
import java.util.Set;


public class EHMM extends ConfigModule {
    public final ExternalConfig<File> bamDirectory = new ExternalConfig<>(File.class);
    public final InternalConfig<Set<String>> antigenClassKeys = new InternalConfig<>(
            Set.of("ATAC-Seq", "Histone_H3K27ac", "Histone_H3K4me1", "Histone_H3K4me3"));
    public final InternalConfig<String> threshold = new InternalConfig<>("5");
    public final InternalConfig<Integer> nStates = new InternalConfig<>(10);
    public final InternalConfig<Double> pseudoCount = new InternalConfig<>(1.0);
    public final InternalConfig<Integer> nBins = new InternalConfig<>(100);
    public final InternalConfig<Integer> nThreads = new InternalConfig<>(21);
    public final InternalConfig<Integer> nSamples = new InternalConfig<>(300);
    public final InternalConfig<Integer> genomicRegionSize = new InternalConfig<>(2000);
}
