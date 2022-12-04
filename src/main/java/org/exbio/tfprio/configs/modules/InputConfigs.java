package org.exbio.tfprio.configs.modules;

import org.exbio.pipejar.configs.ConfigModule;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.ExternalConfig;

import java.io.File;

public class InputConfigs extends ConfigModule {
    public final ExternalConfig<File> chipSeq = new ExternalConfig<>(File.class);
    public final ExternalConfig<File> rnaSeq = new ExternalConfig<>(File.class);
    public final ExternalConfig<File> geneIDs = new ExternalConfig<>(File.class);
    public final ExternalConfig<File> geneAnnotationFile = new ExternalConfig<>(File.class);
}
