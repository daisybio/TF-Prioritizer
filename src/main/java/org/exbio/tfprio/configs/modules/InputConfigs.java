package org.exbio.tfprio.configs.modules;

import org.exbio.pipejar.configs.ClassGetter;
import org.exbio.pipejar.configs.ConfigModule;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.ExternalConfig;
import org.exbio.pipejar.configs.ConfigValidators.PositiveIntegerValidator;
import org.exbio.pipejar.configs.ConfigValidators.StringValidator;

import java.io.File;
import java.util.Map;

public class InputConfigs extends ConfigModule {
    public final ExternalConfig<File> peaks = new ExternalConfig<>(File.class);
    public final ExternalConfig<File> rnaSeq = new ExternalConfig<>(File.class);
    public final ExternalConfig<File> geneIDs = new ExternalConfig<>(File.class);
    public final ExternalConfig<File> geneAnnotationFile = new ExternalConfig<>(File.class);
    public final ExternalConfig<Map<String, String>> sameStages = new ExternalConfig<>(ClassGetter.getMap());
    public final ExternalConfig<String> genome = new ExternalConfig<>(String.class);
    public final ExternalConfig<String> biomartSpecies = new ExternalConfig<>(String.class);
    public final ExternalConfig<String> seqType =
            new ExternalConfig<>(String.class, new StringValidator("chip-seq", "atac-seq", "dnase-seq"));
    public final ExternalConfig<Integer> topTargetGenes =
            new ExternalConfig<>(Integer.class, new PositiveIntegerValidator());
}
