package util.Configs.Modules.FileStructure;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class Igv extends AbstractModule
{
    public final Config<File> root = extend(workingDirectory, "09_IGV_screenshots");
    public final Config<File> ownData = extend(root, "01_own_tf_data");
    public final Config<File> session = extend(root, "A_SESSIONS");
    public final Config<File> chipAtlasData = extend(root, "02_chip_atlas_tf_data");
    public final Config<File> chipAtlasChromosomeWideGenomeWideViews =
            extend(root, "03_chip_atlas_chrWide_genomeWide_views");
    public final Config<File> importantLoci = extend(root, "04_important_loci");

    public final Config<File> igvTopLog2fc = extend(root, "05_top_log2fc");
    public final Config<File> igvTopLog2fc_upregulated = extend(igvTopLog2fc, "01_upregulated");
    public final Config<File> igvTopLog2fc_downregulated = extend(igvTopLog2fc, "02_downregulated");

    public final Config<File> igvDcgTargetGenes = extend(root, "06_dcg_target_genes");

    public Igv(File workingDirectory, File sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
