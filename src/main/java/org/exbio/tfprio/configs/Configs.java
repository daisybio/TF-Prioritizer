package org.exbio.tfprio.configs;

import org.exbio.pipejar.configs.ConfigModuleCollection;
import org.exbio.tfprio.configs.modules.*;

public class Configs extends ConfigModuleCollection {
    public static final InputConfigs inputConfigs = new InputConfigs();
    public static final MixOptions mixOptions = new MixOptions();
    public static final DeSeq2 deSeq2 = new DeSeq2();
    public static final TGene tGene = new TGene();
    public static final TEPIC tepic = new TEPIC();
    public static final DYNAMITE dynamite = new DYNAMITE();
    public static final InternalConfigs internalConfigs = new InternalConfigs();
}
