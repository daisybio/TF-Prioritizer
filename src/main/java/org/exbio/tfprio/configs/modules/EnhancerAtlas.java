package org.exbio.tfprio.configs.modules;

import org.exbio.pipejar.configs.ConfigModule;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.InternalConfig;

import java.util.HashMap;
import java.util.Map;

public class EnhancerAtlas extends ConfigModule {
    public final InternalConfig<String> enhancerBaseURLString = new InternalConfig<>("http://www.enhanceratlas.org/data/download/enhancer/");
    public final InternalConfig<String> enhancerHTMLString = new InternalConfig<>("http://www.enhanceratlas.org/downloadv2.php");
    public final InternalConfig<Map<String, String>> enhancerVersionMap = new InternalConfig<>(new HashMap<>() {{
        put("hs", "hg19");
        put("mm", "mm9");
        put("dr", "danRer10");
        put("dm", "dm3");
        put("ce", "ce10");
        put("rn", "rn5");
        put("sc", "sacCer3");
        put("gg", "galGal4");
        put("ss", "susScr3");
    }});
    public final InternalConfig<Map<String, String>>  genomeToKeyMap = new InternalConfig<>(new HashMap<>() {{
        put("hg38", "hs");
        put("GRCh38", "hs");
        put("hg19", "hs");
        put("GRCh37", "hs");
        put("mm10", "mm");
        put("GRCm38", "mm");
        put("mm9", "mm");
        put("GRCm37", "mm");
        put("zv10", "dr");
        put("danRer10", "dr");
        put("GRCz10", "dr");
        put("zv9", "dr");
        put("GRCz9", "dr");
        put("dm6", "dm");
        put("BDGP6", "dm");
        put("ce10", "ce");
        put("WBcel235", "ce");
        put("rn6", "rn");
        put("Rnor_6.0", "rn");
        put("Galgal4", "gg");
        put("galGal4", "gg");
        put("susScr3", "ss");
        put("Sscrofa10.2", "ss");
        put("sacCer3", "sc");
        put("R64-1-1", "sc");
    }});
}
