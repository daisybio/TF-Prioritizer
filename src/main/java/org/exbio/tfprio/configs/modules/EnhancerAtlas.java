package org.exbio.tfprio.configs.modules;

import org.exbio.pipejar.configs.ConfigModule;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.InternalConfig;

import java.util.Map;

public class EnhancerAtlas extends ConfigModule {
    public final InternalConfig<String> enhancerBaseURLString = new InternalConfig<>("http://www.enhanceratlas.org/data/download/enhancer/");
    public final InternalConfig<String> enhancerHTMLString = new InternalConfig<>("http://www.enhanceratlas.org/downloadv2.php");
    public final InternalConfig<Map<String, String>> enhancerVersionMap = new InternalConfig<>(Map.ofEntries(
            Map.entry("hs", "hg19"),
            Map.entry("mm", "mm9"),
            Map.entry("dr", "danRer10"),
            Map.entry("dm", "dm3"),
            Map.entry("ce", "ce10"),
            Map.entry("rn", "rn5"),
            Map.entry("sc", "sacCer3"),
            Map.entry("gg", "galGal4"),
            Map.entry("ss", "susScr3")
    ));
    public final InternalConfig<Map<String, String>>  genomeToKeyMap = new InternalConfig<>(Map.ofEntries(
            Map.entry("hg38", "hs"),
            Map.entry("GRCh38", "hs"),
            Map.entry("hg19", "hs"),
            Map.entry("GRCh37", "hs"),
            Map.entry("mm10", "mm"),
            Map.entry("GRCm38", "mm"),
            Map.entry("mm9", "mm"),
            Map.entry("GRCm37", "mm"),
            Map.entry("zv10", "dr"),
            Map.entry("danRer10", "dr"),
            Map.entry("GRCz10", "dr"),
            Map.entry("zv9", "dr"),
            Map.entry("GRCz9", "dr"),
            Map.entry("dm6", "dm"),
            Map.entry("BDGP6", "dm"),
            Map.entry("ce10", "ce"),
            Map.entry("WBcel235", "ce"),
            Map.entry("rn6", "rn"),
            Map.entry("Rnor_6.0", "rn"),
            Map.entry("Galgal4", "gg"),
            Map.entry("galGal4", "gg"),
            Map.entry("susScr3", "ss"),
            Map.entry("Sscrofa10.2", "ss"),
            Map.entry("sacCer3", "sc"),
            Map.entry("R64-1-1", "sc")
    ));
}
