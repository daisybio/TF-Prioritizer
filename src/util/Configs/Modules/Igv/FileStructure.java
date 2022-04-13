package util.Configs.Modules.Igv;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class FileStructure extends AbstractModule
{
    public final Config<String> s_session = new Config<>("session.xml");
    public final Config<File> d_root = extend(workingDirectory, "09_IGV_screenshots");
    public final Config<File> d_ownData = extend(d_root, "01_own_tf_data");
    public final Config<File> d_session = extend(d_root, "A_SESSIONS");
    public final Config<File> d_chipAtlasData = extend(d_root, "02_chip_atlas_tf_data");
    public final Config<File> d_chipAtlasChromosomeWideGenomeWideViews =
            extend(d_root, "03_chip_atlas_chrWide_genomeWide_views");
    public final Config<File> d_importantLoci = extend(d_root, "04_important_loci");

    public final Config<File> d_igvTopLog2fc = extend(d_root, "05_top_log2fc");
    public final Config<String> s_igvTopLog2fc_upregulated = new Config<>("01_upregulated");
    public final Config<String> s_igvTopLog2fc_downregulated = new Config<>("02_downregulated");

    public final Config<File> d_igvDcgTargetGenes = extend(d_root, "06_dcg_target_genes");

    public FileStructure(Config<File> workingDirectory, Config<File> sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
