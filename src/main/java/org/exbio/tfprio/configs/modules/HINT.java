package org.exbio.tfprio.configs.modules;

import org.exbio.pipejar.configs.ConfigModule;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.ExternalConfig;

import java.io.File;

public class HINT extends ConfigModule {
    public final ExternalConfig<File> bam_directory = new ExternalConfig<>(File.class);
    public final ExternalConfig<Boolean> paired = new ExternalConfig<>(Boolean.class);
    public final ExternalConfig<String> genome = new ExternalConfig<>(String.class);
}