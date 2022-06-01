package util.Configs.Modules.DeSeq2;

import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.Configs.ConfigTypes.InternalConfig;
import util.Configs.ConfigTypes.SourceDirectoryFileStructure;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.lang.reflect.InvocationTargetException;

import static util.FileManagement.extend;

public class FileStructure extends AbstractModule
{
    public final GeneratedFileStructure d_preprocessing = extend(workingDirectory, "02_A_DESeq2_preprocessing");

    public final GeneratedFileStructure f_mapping = extend(d_preprocessing, "ENSG_SYMBOL_MAP.csv");
    public final GeneratedFileStructure f_mappingScript = extend(d_preprocessing, "ENSG_SYMBOL_MAP.R");
    public final InternalConfig<String> s_preprocessing_MeanCounts = new InternalConfig<>(".tsv");

    public final GeneratedFileStructure d_preprocessing_single = extend(d_preprocessing, "single");
    public final GeneratedFileStructure d_preprocessing_combined = extend(d_preprocessing, "combined");
    public final GeneratedFileStructure d_preprocessing_combinedOriginal = extend(d_preprocessing, "combined_original");
    public final GeneratedFileStructure d_preprocessing_meanCounts =
            extend(d_preprocessing, "symbols_ensg_mean_counts");

    public final GeneratedFileStructure d_preprocessing_tpm = extend(d_preprocessing, "tpm_mapping");
    public final GeneratedFileStructure f_preprocessing_tpm_geneLengths =
            extend(d_preprocessing_tpm, "gene_lengths.tsv");
    public final GeneratedFileStructure d_preprocessing_tpm_scripts = extend(d_preprocessing_tpm, "01_scripts");
    public final GeneratedFileStructure d_preprocessing_tpm_tpmResults = extend(d_preprocessing_tpm, "02_tpm_mappings");
    public final GeneratedFileStructure d_preprocessing_tpm_updated = extend(d_preprocessing_tpm, "03_tpm_updated");

    public final GeneratedFileStructure d_preprocessing_genePositions = extend(d_preprocessing, "gene_positions");
    public final GeneratedFileStructure f_preprocessing_genePositions_script =
            extend(d_preprocessing_genePositions, "get_gene_positions.R");
    public final GeneratedFileStructure f_preprocessing_genePositions_genePositionsPrev =
            extend(d_preprocessing_genePositions, "gene_positions_prev.tsv");
    public final GeneratedFileStructure f_preprocessing_genePositions_version =
            extend(d_preprocessing_genePositions, "version.csv");
    public final GeneratedFileStructure f_preprocessing_genePositions_uplift =
            extend(d_preprocessing_genePositions, "uplift_positions.py");
    public final GeneratedFileStructure f_preprocessing_genePositions_data =
            extend(d_preprocessing_genePositions, "gene_positions.csv");
    public final InternalConfig<String> s_preprocessing_genePositions_mergedEnhancerDbs_bedSuffix =
            new InternalConfig<>(".bed");
    public final InternalConfig<String> s_preprocessing_genePositions_mergedEnhancerDbs_bedFormat =
            new InternalConfig<>("chrom\tchromStart\tchromEnd\tname\treference_genome");
    public final GeneratedFileStructure f_preprocessing_genePositions_mergedEnhancerDbs =
            extend(d_preprocessing_genePositions,
                    "mergedEnhancerDBs" + s_preprocessing_genePositions_mergedEnhancerDbs_bedSuffix.get());


    public final GeneratedFileStructure d_preprocessing_genePositions_enhancerDBs =
            extend(d_preprocessing_genePositions, "enhancerDBs_prev");

    public final GeneratedFileStructure d_outputRaw = extend(workingDirectory, "02_B_DESeq2_output_raw");
    public final GeneratedFileStructure d_rScripts = extend(workingDirectory, "02_C_DESeq2_R_scripts");
    public final GeneratedFileStructure d_output = extend(workingDirectory, "02_D_DESeq2_output");

    public FileStructure(GeneratedFileStructure workingDirectory, SourceDirectoryFileStructure sourceDirectory,
                         Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
