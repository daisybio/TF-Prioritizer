package util.Configs.Modules.Igv;

import util.Configs.ConfigTypes.*;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class FileStructure extends AbstractModule
{
    public final InternalConfig<String> s_session = new InternalConfig<>("session.xml");
    public final GeneratedFileStructure d_root = extend(workingDirectory, "09_IGV_screenshots");
    public final GeneratedFileStructure d_ownData = extend(d_root, "01_own_tf_data");
    public final GeneratedFileStructure d_session = extend(d_root, "A_SESSIONS");
    public final GeneratedFileStructure d_chipAtlasData = extend(d_root, "02_chip_atlas_tf_data");
    public final GeneratedFileStructure d_chipAtlasChromosomeWideGenomeWideViews =
            extend(d_root, "03_chip_atlas_chrWide_genomeWide_views");
    public final GeneratedFileStructure d_importantLoci = extend(d_root, "04_important_loci");

    public final GeneratedFileStructure d_igvTopLog2fc = extend(d_root, "05_top_log2fc");
    public final InternalConfig<String> s_igvTopLog2fc_upregulated = new InternalConfig<>("01_upregulated");
    public final InternalConfig<String> s_igvTopLog2fc_downregulated = new InternalConfig<>("02_downregulated");

    public final GeneratedFileStructure d_igvDcgTargetGenes = extend(d_root, "06_dcg_target_genes");

    public FileStructure(GeneratedFileStructure workingDirectory, SourceDirectoryFileStructure sourceDirectory,
                         Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
