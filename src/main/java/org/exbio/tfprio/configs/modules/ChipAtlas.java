package org.exbio.tfprio.configs.modules;

import org.exbio.pipejar.configs.ClassGetter;
import org.exbio.pipejar.configs.ConfigModule;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.ExternalConfig;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.InternalConfig;
import org.exbio.pipejar.configs.ConfigValidators.ListNotEmptyValidator;

import java.util.List;

public class ChipAtlas extends ConfigModule {
    public final ExternalConfig<Boolean> enabled = new ExternalConfig<>(Boolean.class);
    public final ExternalConfig<List<String>> tissueTypes =
            new ExternalConfig<>(ClassGetter.getList(), new ListNotEmptyValidator<>());
    public final InternalConfig<String> listUrl =
            new InternalConfig<>("http://togodb.biosciencedbc.jp/togodb/release/chip_atlas_file_list.csv");
    public final InternalConfig<String> genomeVersionColName = new InternalConfig<>("genome_assembly");
    public final InternalConfig<String> antigenClassColName = new InternalConfig<>("antigen_class");
    public final InternalConfig<String> antigenRegexColName = new InternalConfig<>("antigen");
    public final InternalConfig<String> cellTypeClassColName = new InternalConfig<>("cell_type_class");
    public final InternalConfig<String> thresholdColName = new InternalConfig<>("threshold");
    public final InternalConfig<String> urlColName = new InternalConfig<>("file_url");
    public final ExternalConfig<String> genomeVersion = new ExternalConfig<>(String.class);
}
