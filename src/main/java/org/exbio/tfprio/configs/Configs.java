package org.exbio.tfprio.configs;

import org.exbio.pipejar.configs.ConfigModuleCollection;
import org.exbio.tfprio.configs.modules.*;

public class Configs extends ConfigModuleCollection {
    public final InputConfigs inputConfigs = new InputConfigs();
    public final MixOptions mixOptions = new MixOptions();
    public final DeSeq2 deSeq2 = new DeSeq2();
    public final TGene tGene = new TGene();
    public final TEPIC tepic = new TEPIC();
    public final DYNAMITE dynamite = new DYNAMITE();
    public final InternalConfigs internalConfigs = new InternalConfigs();
    public final DistributionAnalysis distributionAnalysis = new DistributionAnalysis();
    public final ChipAtlas chipAtlas = new ChipAtlas();
    public final IGV igv = new IGV();
    public final HINT hint = new HINT();
    public final EnhancerAtlas enhancerAtlas = new EnhancerAtlas();
    public final EPDNew epdNew = new EPDNew();
    public final EHMM ehmm = new EHMM();
}
