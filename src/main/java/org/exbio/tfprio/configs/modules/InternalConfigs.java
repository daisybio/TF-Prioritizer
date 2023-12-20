package org.exbio.tfprio.configs.modules;

import org.exbio.pipejar.configs.ConfigModule;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.InternalConfig;

import java.util.HashMap;
import java.util.Map;

public class InternalConfigs extends ConfigModule {
    public final InternalConfig<Map<String, String>> grcSynonymDict = new InternalConfig<>(new HashMap<>() {{
        put("GRCm39", "mm39");
        put("GRCm38", "mm10");
        put("GRCm37", "mm9");
        put("GRCh38", "hg38");
        put("GRCh37", "hg19");
        put("GRCh38.p13", "hg38");
    }});

    public final InternalConfig<String> jasparURL = new InternalConfig<>(
            "https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_non-redundant_pfms_jaspar.txt");

    public final InternalConfig<String> jasparLogoUrl =
            new InternalConfig<>("https://jaspar.elixir.no/static/logos/all/svg/");
}
