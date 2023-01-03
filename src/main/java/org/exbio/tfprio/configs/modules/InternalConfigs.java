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
    }});

    public final InternalConfig<String> jasparURL = new InternalConfig<>(
            "https://jaspar.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_redundant_pfms_jaspar.txt");

    public final InternalConfig<String> jasparLogoUrl =
            new InternalConfig<>("https://jaspar.genereg.net/static/logos/all/svg/");
}
