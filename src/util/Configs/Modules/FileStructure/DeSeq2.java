package util.Configs.Modules.FileStructure;

import util.Configs.Config;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.io.File;
import java.lang.reflect.InvocationTargetException;

public class DeSeq2 extends AbstractModule
{
    public final Config<File> preprocessing =
            new Config<>(new File(workingDirectory.getAbsolutePath() + File.separator + "02_A_DESeq2_preprocessing"));
    public final Config<File> preprocessing_single = extend(preprocessing, "single");
    public final Config<File> preprocessing_combined = extend(preprocessing, "combined");
    public final Config<File> preprocessing_combinedOriginal = extend(preprocessing, "combined_original");
    public final Config<File> preprocessing_geneSymbols = extend(preprocessing, "symbols_ensg_mean_counts");

    public final Config<File> preprocessing_tpm = extend(preprocessing, "tpm_mapping");
    public final Config<File> preprocessing_tpm_scripts = extend(preprocessing_tpm, "01_scripts");
    public final Config<File> preprocessing_tpm_tpmResults = extend(preprocessing_tpm, "02_tpm_mappings");

    public final Config<File> preprocessing_genePositions = extend(preprocessing, "gene_positions");
    public final Config<File> preprocessing_genePositions_enhancerDBs =
            extend(preprocessing_genePositions, "enhancerDBs_prev");

    public final Config<File> outputRaw = extend(workingDirectory, "02_B_DESeq2_output_raw");
    public final Config<File> rScripts = extend(workingDirectory, "02_C_DESeq2_R_scripts");
    public final Config<File> output = extend(workingDirectory, "02_D_DESeq2_output");

    public DeSeq2(File workingDirectory, File sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
