package util.Configs.Modules.Tgene;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class FileStructure extends AbstractModule
{
    public final Config<File> d_root = extend(workingDirectory, "04_TGEN");
    public final Config<File> d_preprocessing = extend(d_root, "01_preprocessing");
    public final Config<File> d_preprocessing_gtf = extend(d_preprocessing, "01_GTF");
    public final Config<String> s_preprocessing_gtf = new Config<>("_transcripts_only.gtf");

    public final Config<File> d_preprocessing_binaryTrees = extend(d_preprocessing, "02_BINARY_TREES");
    public final Config<File> d_preprocessing_binaryTrees_unmerged = extend(d_preprocessing_binaryTrees, "unmerged");
    public final Config<File> d_preprocessing_binaryTrees_merged = extend(d_preprocessing_binaryTrees, "merged");
    public final Config<File> d_preprocessing_binaryTrees_sorted = extend(d_preprocessing_binaryTrees, "sorted");
    public final Config<File> d_output = extend(d_root, "02_output");
    public final Config<String> s_output_links = new Config<>("links.tsv");

    public final Config<File> d_merged = extend(d_root, "03_merged");
    public final Config<File> d_groups = extend(d_root, "04_groups");
    public final Config<String> s_groupts_mergedGroupts = new Config<>("tgene_merged_groups.txt");

    public final Config<File> d_filteredTargetGenes = extend(d_root, "05_filtered_target_genes");
    public final Config<File> d_integrate = extend(d_root, "06_integrate_affinities_self_regulatory");

    public FileStructure(Config<File> workingDirectory, Config<File> sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
