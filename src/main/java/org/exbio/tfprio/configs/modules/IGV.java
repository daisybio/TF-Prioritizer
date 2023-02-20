package org.exbio.tfprio.configs.modules;

import org.exbio.pipejar.configs.ClassGetter;
import org.exbio.pipejar.configs.ConfigModule;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.ExternalConfig;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.InternalConfig;
import org.exbio.pipejar.configs.ConfigValidators.PositiveIntegerValidator;

import java.io.File;
import java.util.List;

public class IGV extends ConfigModule {
    public final InternalConfig<File> igvLib = new InternalConfig<>(new File(System.getenv("IGV"), "lib"));
    public final InternalConfig<File> igvCacheDirectory = new InternalConfig<>(new File(System.getenv("IGV_CACHE")));
    public final InternalConfig<Integer> windowExtend = new InternalConfig<>(50000);
    public final ExternalConfig<File> experimentalFiles = new ExternalConfig<>(File.class);
    public final ExternalConfig<File> signalFiles = new ExternalConfig<>(File.class);
    public final ExternalConfig<List<String>> importantLoci = new ExternalConfig<>(ClassGetter.getList());
    public final ExternalConfig<Integer> topLog2FoldChange =
            new ExternalConfig<>(Integer.class, new PositiveIntegerValidator());
}
