package org.exbio.tfprio.configs.modules;

import org.exbio.pipejar.configs.ConfigModule;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.InternalConfig;

import java.util.Map;

public class EPDNew extends ConfigModule {
    public final InternalConfig<String> baseURL = new InternalConfig<>("ftp://ccg.epfl.ch/epdnew/");
    public final InternalConfig<Map<String, String>> genomeMap = new InternalConfig<>(Map.ofEntries(
            Map.entry("hg38", "H_sapiens/006/Hs_EPDnew_006_hg38.bed"),
            Map.entry("hg19", "H_sapiens/006/Hs_EPDnew_006_hg19.bed"),
            Map.entry("mm10", "M_musculus/003/Mm_EPDnew_003_mm10.bed"),
            Map.entry("mm9", "M_musculus/002/Mm_EPDnew_002_mm9.bed"),
            Map.entry("dr7", "D_rerio/001/Dr_EPDnew_001_danRer7.bed"),
            Map.entry("ce6", "C_elegans/001/Ce_EPDnew_001_ce6.bed"),
            Map.entry("dm6", "D_melanogaster/005/Dm_EPDnew_005_dm6.bed"),
            Map.entry("rn6", "R_norvegicus/001/Rn_EPDnew_001_rn6.bed"),
            Map.entry("sacCer3", "S_cerevisiae/002/Sc_EPDnew_002_sacCer3.bed")
    ));
}
