package util.Configs.Modules.Tgene;

import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.Configs.ConfigTypes.InputFileStructure;
import util.Configs.ConfigTypes.InternalConfig;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class FileStructure extends AbstractModule
{
    public final GeneratedFileStructure d_root = extend(workingDirectory, "04_TGENE");
    public final GeneratedFileStructure d_preprocessing = extend(d_root, "01_preprocessing");
    public final GeneratedFileStructure d_preprocessing_gtf = extend(d_preprocessing, "01_GTF");
    public final GeneratedFileStructure f_transcripts_gtf = extend(d_preprocessing, "transcripts.gtf");
    public final InternalConfig<String> s_preprocessing_gtf = new InternalConfig<>("_transcripts_only.gtf");

    public final GeneratedFileStructure f_preprocessing_regions = extend(d_preprocessing, "regions.tsv");
    public final GeneratedFileStructure d_output = extend(d_root, "02_output");
    public final InternalConfig<String> s_output_links = new InternalConfig<>("links.tsv");

    public final GeneratedFileStructure d_merged = extend(d_root, "03_merged");
    public final GeneratedFileStructure d_groups = extend(d_root, "04_groups");
    public final InternalConfig<String> s_groups_mergedGroups = new InternalConfig<>("tgene_merged_groups.txt");

    public final GeneratedFileStructure d_filteredTargetGenes = extend(d_root, "05_filtered_target_genes");
    public final GeneratedFileStructure d_integrate = extend(d_root, "06_integrate_affinities_self_regulatory");

    public FileStructure(GeneratedFileStructure workingDirectory, InputFileStructure sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
