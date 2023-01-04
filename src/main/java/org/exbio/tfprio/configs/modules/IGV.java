package org.exbio.tfprio.configs.modules;

import org.exbio.pipejar.configs.ConfigModule;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.ExternalConfig;

import java.io.File;

public class IGV extends ConfigModule {
    public final ExternalConfig<File> igvTools = new ExternalConfig<>(File.class);
}
