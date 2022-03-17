package util.Configs.Modules.DeSeq2;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class FileStructure extends AbstractModule
{
    public final Config<File> d_preprocessing = extend(workingDirectory, "02_A_DESeq2_preprocessing");

    public final Config<File> f_mapping = extend(d_preprocessing, "ENSG_SYMBOL_MAP.csv");
    public final Config<File> f_mappingScript = extend(d_preprocessing, "ENSG_SYMBOL_MAP.R");
    public final Config<String> s_preprocessing_MeanCounts = new Config<>("_meanCounts.txt");

    public final Config<File> d_preprocessing_single = extend(d_preprocessing, "single");
    public final Config<File> d_preprocessing_combined = extend(d_preprocessing, "combined");
    public final Config<File> d_preprocessing_combinedOriginal = extend(d_preprocessing, "combined_original");
    public final Config<File> d_preprocessing_meanCounts = extend(d_preprocessing, "symbols_ensg_mean_counts");

    public final Config<File> d_preprocessing_tpm = extend(d_preprocessing, "tpm_mapping");
    public final Config<File> f_preprocessing_tpm_geneLengths = extend(d_preprocessing_tpm, "gene_lengths.tsv");
    public final Config<File> d_preprocessing_tpm_scripts = extend(d_preprocessing_tpm, "01_scripts");
    public final Config<String> s_preprocessing_tpm_scripts = new Config<>("_get_tpms.R");
    public final Config<File> d_preprocessing_tpm_tpmResults = extend(d_preprocessing_tpm, "02_tpm_mappings");
    public final Config<String> s_preprocessing_tpm_mappings = new Config<>("_tpms.csv");

    public final Config<File> d_preprocessing_genePositions = extend(d_preprocessing, "gene_positions");
    public final Config<File> f_preprocessing_genePositions_script =
            extend(d_preprocessing_genePositions, "get_gene_positions.R");
    public final Config<File> f_preprocessing_genePositions_genePositionsPrev =
            extend(d_preprocessing_genePositions, "gene_positions_prev.csv");
    public final Config<File> f_preprocessing_genePositions_version =
            extend(d_preprocessing_genePositions, "version.csv");
    public final Config<File> f_preprocessing_genePositions_uplift =
            extend(d_preprocessing_genePositions, "uplift_positions.py");
    public final Config<File> f_preprocessing_genePositions_data =
            extend(d_preprocessing_genePositions, "gene_positions.csv");
    public final Config<String> s_preprocessing_genePositions_mergedEnhancerDbs_bedSuffix = new Config<>(".bed", true);
    public final Config<String> s_preprocessing_genePositions_mergedEnhancerDbs_bedFormat =
            new Config<>("chrom\tchromStart\tchromEnd\tname\treference_genome", true);
    public final Config<File> f_preprocessing_genePositions_mergedEnhancerDbs = extend(d_preprocessing_genePositions,
            "mergedEnhancerDBs" + s_preprocessing_genePositions_mergedEnhancerDbs_bedSuffix.get());


    public final Config<File> d_preprocessing_genePositions_enhancerDBs =
            extend(d_preprocessing_genePositions, "enhancerDBs_prev");

    public final Config<File> d_outputRaw = extend(workingDirectory, "02_B_DESeq2_output_raw");
    public final Config<File> d_rScripts = extend(workingDirectory, "02_C_DESeq2_R_scripts");
    public final Config<File> d_output = extend(workingDirectory, "02_D_DESeq2_output");
    public final Config<String> s_output_dynamite = new Config<>("_DYNAMITE.tsv");

    public FileStructure(Config<File> workingDirectory, Config<File> sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
