package util;

import java.io.File;

public class Options_intern
{
    /**
     * COM2POSE OPTIONS intern
     */
    public String config_data_path="";
    public String com2pose_working_directory="";
    public String path_to_COM2POSE = "";

    public boolean write_to_logfile = true;

    /**
     * COM2POSE private options - cannot be set from the outside
     */
    public String folder_name_deseq2_preprocessing = "DESeq2_preprocessing";
    public String folder_name_deseq2_preprocessing_single = "single";
    public String folder_name_deseq2_preprocessing_combined = "combined";
    public String folder_name_deseq2_output_raw = "DESeq2_output_raw";
    public String folder_name_deseq2_output = "DESeq2_output";
    public String folder_name_deseq2_R_scripts = "DESeq2_R_scripts";
    public String folder_name_tepic_output_raw = "TEPIC_output_raw";
    public String folder_name_tepic_postprocessing = "TEPIC_postprocessing";
    public String folder_name_tepic_postprocessing_input = "input";
    public String folder_name_tepic_postprocessing_output = "output";
    public String folder_name_tepic_postprocessing_output_mean_affinities= "MeanAffinities";
    public String folder_name_tepic_postprocessing_output_ratios = "Ratios";
    public String folder_output_preprocessing_DYNAMITE = "DYNAMITE_preprocessing";
    public String folder_output_preprocessing_DYNAMITE_integrateData = "integrateData";
    public String folder_output_preprocessing_DYNAMITE_prepareClass = "prepareClassification";
    public String folder_out_put_DYNAMITE = "DYNAMITE_output";


    public String file_suffix_deseq2_preprocessing_meanCounts = "_meanCounts.txt";
    public String file_suffix_deseq2_output_DYNAMITE = "_DYNAMITE.tsv";
    public String file_suffix_tepic_postprocessing_output_mean_affinities ="Mean_Affinities_";
    public String file_suffix_tepic_postprocessing_output_ratios ="Ratio_Affinities_";
    public String file_suffix_output_preprocessing_DYNAMITE_integrateData_log2coeff= "Integrated_Data_Log2_Quotient.txt";
    public String file_suffix_output_preprocessing_DYNAMITE_prepClass ="Integrated_Data_For_Classification.txt";


    public String directory_for_tepic_scripts = "ext"+ File.separator+"TEPIC"+File.separator+"TEPIC";
    public String directory_for_tepic_scripts_code = directory_for_tepic_scripts+File.separator+"Code";
    public String directory_for_tepic_scripts_code_tepic_sh = directory_for_tepic_scripts_code + File.separator+"TEPIC.sh";
    public String directory_for_tepic_DYNAMITE = directory_for_tepic_scripts+File.separator+"MachineLearningPipelines"+File.separator+"DYNAMITE"+File.separator+"Scripts";


    /*
    ######################
    ##DESeq2 parameters###
    ######################*/
    //#[REQ]: count files in nfcore RNA-seq format (each line is a count of one gene), directory must be ordered like: TP1 - samples_1,...,samples_n;....
    public String deseq2_input_directory="";
    //#[REQ]: gene ID file from nfcore RNA-seq (each line names one gene ID (ENSG) - must be same order as deseq2_input_directory files
    public String deseq2_input_gene_id="";
    //#[OPT]: minimum count over all samples of two timepoints for DESeq2, default: 0
    public int deseq2_count_threshold=0;
    /*
    ######################
    ##TEPIC parameters####
    ######################*/
    //#[REQ]: PEAK files (like from MACS), directory must be orderd like: Timepoint1 - HistoneModification1 - samples_1,...,sample_n;...;...
    public String tepic_input_directory="";
    //#[REQ]:fasta files (in RefSeq format, without \"chr\" prefix)
    public String tepic_input_ref_genome="";
    //#[REQ]: path to position specific energy matrix used for TRAP (different matrices can be found in ~/COM2POSE/ext/TEPIC/TEPIC/PWMs/2.1)
    public String tepic_path_pwms="";
    //#[OPT]: number of cores
    public int tepic_cores=1;
    //#[OPT]: bedgraph file containing open chromatin signal, e.g. DNase1-seq, or Histone-Mark signal
    public String tepic_bed_chr_sign="";
    //#[OPT]: column in the tepic_input_directory file containg the average per base signal within a peak. If this option is used, the tepic_bed_chr_sign option must not be used
    public int tepic_column_bedfile=-1;
    //#[OPT]: gene annotation file, required to generate the gene view, required for TPM filter
    public String tepic_gene_annot="";
    //#[OPT]: size of the window to be considered to generate gene view (default 50kb)]
    public int tepic_window_size=50000;
    //#[OPT]: path to annotation file and annotate only DNase peaks that are within a window specified by the tepic_window_size option around all genes contained in the gene annotation file specified by this option
    public String tepic_onlyDNasePeaks="";
    //#[OPT]: indicating whether exponential decay should be used (default TRUE)
    public boolean tepic_exponential_decay=true;
    //#[OPT]: flag to be set if affinities should not be normalised by peak length
    public boolean tepic_not_norm_peak_length=false;
    //#[OPT]: flag to be set if peak features for peak length and peak counts should not be generated
    public boolean tepic_not_generated=false;
    //#[OPT]: if tepic_bed_chr_sign or tepic_column_bedfile is used together with this flag, the original (Decay-)Scaling formulation of TEPIC is used to compute gene-TF scores
    public boolean tepic_original_decay=false;
    //#[OPT]: path to a tab delimited file containing the length of the used PSEMs
    public String tepic_psems_length="";
    //#[OPT]: flag to be set if the entire gene body should be screened for TF binding. The search window is extended by a region half of the size that is specified by the tepic_window_size option upstream of the genes 5' TSS
    public boolean tepic_entire_gene_body=false;
    //#[OPT]: flag indicating that the output of TEPIC should be zipped
    public boolean tepic_zipped=false;
    //#[OPT]: path to a set of background sequences that should be used to compute to generate a binary score for TF binding. Mutually exclusive to the tepic_2bit option
    public String tepic_background_seq="";
    //#[OPT]: path to a 2bit representation of the reference genome, required to generate a binary score for TF binding. The binary score is generated in addition to the standard affinity values. Mutually exclusive to the tepic_background_seq option
    public String tepic_2bit="";
    //#[OPT]: p-value cut off used to determine a cut off to derive a binary score for TF binding (default 0.05)
    public double tepic_pvalue=0.05;
    //#[OPT]: minutes that should be spend at most per chromosome to find matching random regions (default 3)
    public int tepic_minutes_per_chr=3;
    //#[OPT]: flag indicating that the reference genome contains a chr prefix
    public boolean tepic_chr_prefix=false;
    //#[OPT]: flag indicating that the annotation should be transcript and not gene based
    public boolean tepic_transcript_based=false;
    //#[OPT]: loop list file containing chromatin contacts
    public String tepic_loop_list="";
    //#[OPT]: size of the loop window used around a genes promoter to link chromatin loops to genes (default 5000)
    public int tepic_loop_windows=5000;
    //#[OPT]: parameter to be set if only peak features should be computed (default FALSE)
    public boolean tepic_only_peak_features=false;
    //#[OPT]: set (T)ranscripts (P)er (M)illion cutoff, default: TPM filter not active
    public double tepic_tpm_cutoff=-1;
    //#[OPT]: path to input file ensg to gene symbol file, required for TPM filter
    public String tepic_ensg_symbol="";

   /*######################
     ##DYNAMITE parameters#
     ######################*/
   //####preprocessing####
    //#[OPT] ##DO NOT SET IF FULL PIPELINE IS RUN: Position of the gene IDs in the expression file",default=0
    public int dynamite_preprocessing_integrate_data_geneIDs=0;
    //#[OPT] ##DO NOT SET IF FULL PIPELINE IS RUN: Position of the gene expression estimate in the expression file, default=1"
    public int dynamite_preprocessing_integrate_data_log2fc=1;
    //#[OPT] File containing gene IDs that should be considered
    public String dynamite_preprocessing_integrate_data_consider_geneFile="";
    //####DYNAMITE####
    // #[REQ]: Name of the response variable (e.g. Expression)
    public String dynamite_out_var="";
    //#[OPT]: Number of the cores to use (1 as default)
    public int dynamite_cores=1;
    //#[OPT]: Alpha parameter stepsize (0.1 as default)
    public double dynamite_alpha=0.1;
    //#[OPT]:Size of test data (0.2 as default)
    public double dynamite_testsize=0.2;
    //#[OPT]: Number of outer folds for model validation (3 as default)
    public int dynamite_Ofolds=3;
    //#[OPT]: Number of inner cross validation folds (6 as default)
    public int dynamite_Ifolds=6;
    //#[OPT]: Flag indicating whether the data should be balanced through downsampling (default TRUE)
    public boolean dynamite_balanced=true;
    //#[OPT]: Flag indicating whether performance measures should be computed (default TRUE)
    public boolean dynamite_performance=true;
    //#[OPT]: Flag indicating whether a model should be learned on randomised data (default FALSE)
    public boolean dynamite_randomise=false;

}
