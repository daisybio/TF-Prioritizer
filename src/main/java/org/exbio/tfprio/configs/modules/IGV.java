package org.exbio.tfprio.configs.modules;

import org.exbio.pipejar.configs.ConfigModule;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.ExternalConfig;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.InternalConfig;

import java.io.File;

public class IGV extends ConfigModule {
    public final InternalConfig<File> igvExecutable = new InternalConfig<>(new File(System.getenv("IGV"), "igv.sh"));
    public final InternalConfig<Integer> windowExtend = new InternalConfig<>(50000);
    public final ExternalConfig<File> experimentalFiles = new ExternalConfig<>(File.class);
    public final ExternalConfig<File> signalFiles = new ExternalConfig<>(File.class);
    public final ExternalConfig<String> genome = new ExternalConfig<>(String.class);
}
