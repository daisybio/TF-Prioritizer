package util.Configs.Modules.FileStructure;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class Tgen extends AbstractModule
{
    public final Config<File> root = extend(workingDirectory, "04_TGEN");
    public final Config<File> preprocessing = extend(root, "01_preprocessing");
    public final Config<File> preprocessing_gtf = extend(preprocessing, "01_GTF");
    public final Config<File> preprocessing_binaryTrees = extend(preprocessing, "02_BINARY_TREES");
    public final Config<File> preprocessing_binaryTrees_unmerged = extend(preprocessing_binaryTrees, "unmerged");
    public final Config<File> preprocessing_binaryTrees_merged = extend(preprocessing_binaryTrees, "merged");
    public final Config<File> preprocessing_binaryTrees_sorted = extend(preprocessing_binaryTrees, "sorted");
    public final Config<File> output = extend(root, "02_output");
    public final Config<File> merged = extend(root, "03_merged");
    public final Config<File> groups = extend(root, "04_groups");
    public final Config<File> filteredTargetGenes = extend(root, "05_filtered_target_genes");
    public final Config<File> integrate = extend(root, "06_integrate_affinities_self_regulatory");

    public Tgen(File workingDirectory, File sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
