package util.Configs.Modules.DeSeq2;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class FileStructure extends AbstractModule
{
    public final Config<File> d_preprocessing =
            new Config<>(new File(workingDirectory.getAbsolutePath() + File.separator + "02_A_DESeq2_preprocessing"));
    public final Config<File> d_preprocessing_single = extend(d_preprocessing, "single");
    public final Config<File> d_preprocessing_combined = extend(d_preprocessing, "combined");
    public final Config<File> d_preprocessing_combinedOriginal = extend(d_preprocessing, "combined_original");
    public final Config<File> d_preprocessing_geneSymbols = extend(d_preprocessing, "symbols_ensg_mean_counts");

    public final Config<File> d_preprocessing_tpm = extend(d_preprocessing, "tpm_mapping");
    public final Config<File> d_preprocessing_tpm_scripts = extend(d_preprocessing_tpm, "01_scripts");
    public final Config<File> d_preprocessing_tpm_tpmResults = extend(d_preprocessing_tpm, "02_tpm_mappings");

    public final Config<File> d_preprocessing_genePositions = extend(d_preprocessing, "gene_positions");
    public final Config<File> d_preprocessing_genePositions_enhancerDBs =
            extend(d_preprocessing_genePositions, "enhancerDBs_prev");

    public final Config<File> d_outputRaw = extend(workingDirectory, "02_B_DESeq2_output_raw");
    public final Config<File> d_rScripts = extend(workingDirectory, "02_C_DESeq2_R_scripts");
    public final Config<File> d_output = extend(workingDirectory, "02_D_DESeq2_output");

    public FileStructure(File workingDirectory, File sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
