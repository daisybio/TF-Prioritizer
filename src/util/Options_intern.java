package util;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

public class Options_intern
{
    /**
     * COM2POSE OPTIONS intern
     */
    public String config_data_path="";
    public String com2pose_working_directory="";
    public String path_to_COM2POSE = "";

    public String path_tgen = "";

    public boolean write_to_logfile = true;
    public boolean do_ensg_mapping= true;
    public boolean calculate_tpm_lengths = true;
    public boolean calculcate_gene_positions=true;

    /**
     * COM2POSE private options - cannot be set from the outside
     */

    public String folder_name_usual_working_dir_name = "working_dir";

    public String distribution_analysis_all_name = "ALL";

    public String folder_name_mix_option = "00_MIX_OPTION";
    public String folder_name_mix_option_sample_mix_preprocessing = "00_PREPROCESSING_SAMPLE_MIX";
    public String folder_name_mix_option_sample_mix = "01_SAMPLE_MIX";
    public String folder_name_mix_option_preprocess_hm_mix = "02_PREPROCESS_HM_MIX";
    public String folder_name_mix_option_hm_mix = "03_HM_MIX";
    public String folder_name_mix_option_mutually_exclusive = "02_MUTUALLY_EXCLUSIVE";
    public String folder_name_mix_option_mutually_exclusive_preprocessing = "01_PREPROCESSING";
    public String folder_name_mix_options_mutually_exclusive_input = "02_NEW_INPUT";
    public String folder_name_blacklisted_regions = "01_blacklisted_regions";
    public String folder_name_blacklisted_regions_preprocessing = "01_preprocessing";
    public String folder_name_blacklisted_regions_preprocessing_perChr = "01_perChr";
    public String folder_name_blacklisted_regions_preprocessing_sorted = "02_sorted";
    public String folder_name_blacklisted_regions_new_input = "02_new_input";
    public String folder_name_deseq2_preprocessing = "02_A_DESeq2_preprocessing";
    public String folder_name_deseq2_preprocessing_single = "single";
    public String folder_name_deseq2_preprocessing_combined = "combined";
    public String folder_name_deseq2_preprocessing_combined_original = "combined_original";
    public String folder_name_deseq2_preprocessing_gene_symbols = "symbols_ensg_mean_counts";
    public String folder_name_deseq2_preprocessing_tpm = "tpm_mapping";
    public String folder_name_deseq2_preprocessing_gene_positions = "gene_positions";
    public String folder_name_deseq2_preprocessing_tpm_scripts = "01_scripts";
    public String folder_name_deseq2_preprocessing_tpm_results = "02_tpm_mappings";
    public String folder_name_deseq2_output_raw = "02_B_DESeq2_output_raw";
    public String folder_name_deseq2_R_scripts = "02_C_DESeq2_R_scripts";
    public String folder_name_deseq2_output = "02_D_DESeq2_output";
    public String folder_name_tepic_output_raw = "03_A0_TEPIC_output_raw";
    public String folder_name_tepic_output_raw_shuffle = "03_A1_TEPIC_output_raw_shuffle";
    public String folder_name_tepic_postprocessing = "03_B_TEPIC_postprocessing";
    public String folder_name_tepic_postprocessing_input = "input";
    public String folder_name_tepic_postprocessing_output = "output";
    public String folder_name_tepic_postprocessing_tfs = "TFs";
    public String folder_name_tepic_postprocessing_output_mean_affinities= "MeanAffinities";
    public String folder_name_tepic_postprocessing_output_ratios = "Ratios";
    public String folder_name_tgen = "04_TGENE";
    public String folder_name_tgen_preprocessing = "01_preprocessing";
    public String folder_name_tgen_preprocessing_gtf = "01_GTF";
    public String folder_name_tgen_preprocessing_binary_trees = "02_BINARY_TREES";
    public String folder_name_tgen_preprocessing_binary_trees_unmerged = "unmerged";
    public String folder_name_tgen_preprocessing_binary_trees_merged = "merged";
    public String folder_name_tgen_preprocessing_binary_trees_sorted = "sorted";
    public String folder_name_tgen_output = "02_output";
    public String folder_name_tgen_merged = "03_merged";
    public String folder_name_tgen_groups = "04_groups";
    public String folder_name_tgen_filter_target_genes = "05_filtered_target_genes";
    public String folder_name_tgen_integrate = "06_integrate_affinities_self_regulatory";
    public String folder_output_preprocessing_DYNAMITE = "05_A_DYNAMITE_preprocessing";
    public String folder_output_preprocessing_DYNAMITE_integrateData = "integrateData";
    public String folder_output_preprocessing_DYNAMITE_prepareClass = "prepareClassification";
    public String folder_output_preprocessing_DYNAMITE_install_required_packages="X_install_required_packages";
    public String folder_out_put_DYNAMITE = "05_B_DYNAMITE_output";
    public String folder_out_data_plots = "06_A_PLOTS_DATA";
    public String folder_plots = "06_B_PLOTS_output";
    public String folder_out_analysis_data = "06_C_PLOTS_ANALYSIS";
    public String folder_out_analysis_data_TP_LEVEL = "01_TP_LEVEL";
    public String folder_out_analysis_data_HM_LEVEL = "02_HM_LEVEL";
    public String folder_out_analysis_data_WEBSITE_OVERVIEW = "03_WEBSITE_OVERVIEW";
    public String folder_out_analysis_data_WEBSITE_NUMBER_TFs = "04_WEBSITE_NUMBER_TFs";
    public String folder_out_target_genes ="06_D_PLOTS_TARGET_GENES";
    public String folder_out_target_genes_all_different = "all_data_different";
    public String folder_out_target_genes_same = "all_data_same";
    public String folder_out_target_genes_dcg = "06_E_TARGET_GENES_DCG";
    public String folder_out_distribution ="07_DISTRIBUTION_ANALYSIS";
    public String folder_out_distribution_analyzed_tfs ="01_ANALYSED_TFS";
    public String folder_out_distribution_tf_tg_scores = "02_TF_TG_SCORES";
    public String folder_out_distribution_tf_tg_scores_background_distr = "01_BACKGROUND_DISTRIBUTION";
    public String folder_out_distribution_tf_tg_scores_background_distr_ALL = "01_ALL";
    public String folder_out_distribution_tf_tg_scores_background_distr_HM = "02_HM";
    public String folder_out_distribution_tf_tg_scores_tf_distributions = "02_TF_DISTRIBUTION";
    public String folder_out_distribution_tf_tg_scores_tf_distributions_ALL = "01_ALL";
    public String folder_out_distribution_tf_tg_scores_tf_distributions_HM = "02_HM";
    public String folder_out_distribution_plots_scripts = "03_PLOT_SCRIPTS";
    public String folder_out_distribution_plots_script_ALL = "01_ALL";
    public String folder_out_distribution_plots_scripts_HM = "02_HM";
    public String folder_out_distribution_plots = "04_DISTR_PLOTS";
    public String folder_out_distribution_plots_ALL = "01_ALL";
    public String folder_out_distribution_plots_HM = "02_HM";
    public String folder_out_distribution_stats = "05_STATS";
    public String folder_out_distribution_stats_ALL = "01_ALL";
    public String folder_out_distribution_stats_HM = "01_HM";
    public String folder_out_distribution_hypergeometric_test = "06_HYPERGEOMETRIC_TEST";
    public String folder_out_distribution_mwu_scripts = "07_SCRIPTS_MANN_WHITNEYU_PLOTS";
    public String folder_out_distribution_mwu_plots = "08_PLOTS_MANN_WHITNEYU_PLOTS";
    public String folder_out_distribution_dcg = "09_DISCOUNTED_CUMULATIVE_GAIN";
    public String folder_out_website="Z_WEBSITE_OVERVIEW";
    public String folder_out_website_interactive_plots = "INTERACTIVE_PLOTS";
    public String folder_out_website_basics = "WEBSITE_BASICS";
    public String folder_out_website_basics_website = "WEBSITE";
    public String folder_out_website_basics_website_css = "CSS";
    public String folder_out_website_basics_website_images = "images";
    public String folder_out_website_htmls_regression_coefficients = "HTMLS_REGRESSION_COEFFICIENT_ANALYSIS";
    public String folder_out_website_htmls_distribution_analysis = "HTMLS_DISTRIBUTION_ANALYSIS";
    public String folder_out_website_htmls_distribution_analysis_ALL = "01_ALL";
    public String folder_out_website_htmls_distribution_analysis_HM = "02_HM";
    public String folder_out_website_htmls_distribution_analysis_cumulative_gain = "03_CUMULATIVE_GAIN";
    public String folder_out_website_htmls_TFs = "TFs";
    public String folder_out_website_interactive_plots_tps="TIMEPOINTS_CONDITIONS";
    public String folder_out_website_interactive_plots_overview="OVERVIEW";
    public String folder_out_website_plots_distribution_analysis = "DISTRIBUTION_ANALYSIS_PLOTS";
    public String folder_out_chip_atlas = "08_CHIP_ATLAS_EVALUATION_PEAK_DATA";
    public String folder_out_chip_atlas_list = "01_CHIP_ATLAS_LIST";
    public String folder_out_chip_atlas_peak_files = "02_CHIP_ATLAS_PEAK_FILES";
    public String folder_out_igv = "09_IGV_screenshots";
    public String folder_out_igv_own_data="01_own_tf_data";
    public String folder_out_igv_chip_atlas_data="02_chip_atlas_tf_data";
    public String folder_out_igv_chipAtlas_chrWide_genomeWide_views = "03_chip_atlas_chrWide_genomeWide_views";


    public String file_suffix_deseq2_mapping = "ENSG_SYMBOL_MAP.csv";
    public String file_suffix_deseq2_preprocessing_meanCounts = "_meanCounts.txt";
    public String file_suffix_deseq2_preprocessing_tpm_mapping_get_geneLengths_script ="get_gene_lengths.R";
    public String file_suffix_deseq2_preprocessing_tpm_mapping_geneLengths_file = "gene_lengths.csv";
    public String file_suffix_deseq2_preprocessing_tpm_mapping_get_tpm_mappings_script = "_get_tpms.R";
    public String file_suffix_deseq2_preprocessing_tpm_mapping_get_tpm_mappings_data = "_tpms.csv";
    public String file_suffix_deseq2_preprocessing_gene_positions_script="get_gene_positions.R";
    public String file_suffix_deseq2_preprocessing_gene_positions_data = "gene_positions.csv";
    public String file_suffix_deseq2_output_DYNAMITE = "_DYNAMITE.tsv";
    public String file_suffix_tepic_output_regions_to_target_genes = "regions_to_target_genes.csv";
    public String file_suffix_tepic_postprocessing_all_tfs = "ALL_TFs.csv";
    public String file_suffix_tepic_postprocessing_output_mean_affinities ="Mean_Affinities_";
    public String file_suffix_tepic_postprocessing_output_ratios ="Ratio_Affinities_";
    public String file_suffix_tepic_postprocessing_tfs_tfs = "tfs.csv";
    public String file_suffix_tgen_preprocess_gtf = "_transcripts_only.gtf";
    public String file_suffix_tgen_output = "links.tsv";
    public String file_suffic_tgen_output_groups = "tgene_merged_groups.txt";
    public String file_suffix_output_preprocessing_DYNAMITE_integrateData_log2coeff= "Integrated_Data_Log2_Quotient.txt";
    public String file_suffix_output_preprocessing_DYNAMITE_prepClass ="Integrated_Data_For_Classification.txt";
    public String file_suffix_output_preprocessing_DYNAMITE_install_required_packages_dynamite="install_required_packages_DYNAMITE.R";
    public String file_suffix_dynamite_output_to_be_plotted = "Regression_Coefficients_Entire_Data_Set_Integrated_Data_For_Classification.txt";
    public String file_suffix_analysis_plot_data_hm_level_different ="all_data_different.csv";
    public String file_suffix_analysis_plot_data_hm_level_same ="all_data_same.csv";
    public String file_suffix_website_analysis_tf_available = "available_tfs.csv";
    public String file_suffix_distribution_analysis_analysed_tfs = "analysable_tfs.csv";
    public String file_suffix_distribution_analysis_distributions = "distribution.csv";
    public String file_suffix_distribution_analysis_python_script = "analyse_tf_to_background.py";
    public String file_suffix_distribution_analysis_plot_stats = "stats.csv";
    public String file_suffix_distribution_analysis_hypergeometric_test_rscript = "HYPERGEOMETRIC_TEST.R";
    public String file_suffix_distribution_analysis_hypergeometric_test_output = "STATS.csv";
    public String file_suffix_distribution_analysis_mann_whitneyU_plot_scripts = "mann_whitneyU_plots.R";
    public String file_suffix_distribution_analysis_dcg = "dcg_stats.csv";
    public String file_suffix_chip_atlas_list_zipped = "chip_atlas_file_list.zip";
    public String file_suffix_chip_atlas_list_csv = "chip_atlas_file_list.csv";

    public String html_report_home_home = "HOME.html";
    public String html_report_home_regression_coefficient_analysis = "HOME_REGRESSION.html";
    public String html_report_home_distribution_analysis = "HOME_DISTRIBUTION.html";
    public String html_report_home_regression_distribution_analysis_all = "ALL.html";

    public String html_report_levels_home = "HOME";
    public String html_report_levels_2_steps= "LEVEL_2";
    public String html_report_levels_3_steps = "THRESHOLD";
    public String html_report_levels_4_steps = "TF";

    public String analysis_types_distribution_analysis ="DISTRIBUTION_ANALYSIS";
    public String analysis_types_regression_coefficient_analysis = "REGRESSION_COEFFICIENT";


    public String directory_for_tepic_scripts = "ext"+ File.separator+"TEPIC"+File.separator+"TEPIC";
    public String directory_for_tepic_scripts_code = directory_for_tepic_scripts+File.separator+"Code";
    public String directory_for_tepic_scripts_code_tepic_sh = directory_for_tepic_scripts_code + File.separator+"TEPIC.sh";
    public String directory_for_tepic_DYNAMITE = directory_for_tepic_scripts+File.separator+"MachineLearningPipelines"+File.separator+"DYNAMITE"+File.separator+"Scripts";

    /*#################################
      ##PREPROCESSING MIX OPTIONS######
      #################################*/
    //#[OPT]: if set a mix of either the Histone Modification Level (HM_LEVEL) or the Sample Level (SAMPLE_LEVEL) will be performed, mix_option is required
    public String mix_level="";
    //#[OPT]: set histone marks or samples will be mixed with option: UNION (all peaks of all HMs / samples will be used), INTERSECTION (only peaks, which are in at least 2 HMs / samples will be used)
    public String mix_option="";
    //#[OPT]: minimal occurence of peaks in intersection, default 2
    public int mix_occurence_intersection=2;
    //#[OPT]: mutually exclusive peaks. Use only peaks, which are mutually exclusive, default: FALSE
    public boolean mix_mutually_exclusive=false;
    //#[OPT]: includes peaks, which are not mutually exclusive, but are differential in peak signal (can only be used if mix_mutually_exclusive=TRUE), default = TRUE
    public boolean mix_mutually_exclusive_diff_peak_signals=true;

    /*
    #########################
    ##blacklist parameters###
    #########################*/
    //#[OPT]: if blacklist filter should be used provide path to BED file here (BED files can be found at https://github.com/Boyle-Lab/Blacklist/tree/master/lists
    public String black_list_dir="";
    //#[OPT]: which regions should be ignored: Low_Mappability;High_Signal_Region
    public HashSet<String> black_list_signals= new HashSet<>(Arrays.asList(new String[]{"Low_Mappability".toUpperCase(), "High_Signal_Region".toUpperCase()}));

    /*
    ######################
    ##DESeq2 parameters###
    ######################*/
    //#[REQ]: count files in nfcore RNA-seq format (each line is a count of one gene), directory must be ordered like: TP1 - samples_1,...,samples_n;....
    public String deseq2_input_directory="";
    //#[REQ]: gene ID file from nfcore RNA-seq (each line names one gene ID (ENSG) - must be same order as deseq2_input_directory files
    public String deseq2_input_gene_id="";
    //#[OPT]:
    public String deseq2_biomart_dataset_species="";
    //#[REQ]: biomart dataset column for gene symbols (e.g. human: hgnc_symbol, mouse: mgi_symbol)
    public String deseq2_biomart_dataset_symbol_column="mgi_symbol";
    //#[OPT]: minimum count over all samples of two timepoints for DESeq2, default: 0
    public int deseq2_count_threshold=0;
    //#[OPT]: TPM filter for RNA-seq data, default: 0.0
    public double deseq2_tpm_filter=0.0;
    /*
    ######################
    ##TEPIC parameters####
    ######################*/
    //#[REQ]: PEAK files (like from MACS), directory must be orderd like: Timepoint1 - HistoneModification1 - samples_1,...,sample_n;...;...
    public String tepic_input_directory="";
    //if black listed regions filter is used the former data link is here
    public String tepic_input_prev ="";
    //original tepic input directory
    public String tepic_input_original = "";
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
    //#[OPT]: use only tgene linked target genes default: true
    public boolean tepic_tgene_target_genes=true;
    //#[OPT]: randomize TEPIC output, default: false
    public boolean tepic_randomize_tf_gene_matrix=false;

    /*######################
     ####TGENE parameters###
     #######################*/
    //#[OPT]:  if no locus is found within window size, the nearest locus is used, default:false
    public boolean tgen_no_closest_locus=false;
    //#[OPT]:  if no tss is found within window size, the nearest locus is used, default:false
    public boolean tgen_no_closest_tss=false;
    //#[OPT]:  window size of tgen, default: 500000
    public int tgen_max_link_distances=500000;
    //#[OPT]: max pvalue, which is accepted, default 0.05
    public double tgen_pvalue=0.05;
    //#[REQ]: should self regulatory TFs {OPT-SELF-REG} be increased? default: false
    public boolean tgen_self_regulatory=false;
    //#[OPT]: {OPT-SELF-REG} value for consus approach (e.g. 0.5 = 50:50 TGene:TEPIC, 0.4 = 40:60 TGene:TEPIC), default: 0.5
    public double tgen_consensus=0.5;
    //#[OPT]: {OPT-SELF-REG} can be "INCREASE_TGENE_TFS" (increases TEPIC TF affinities for TFs found in TGENE at target gene) or "DECREASE_NOT_TGENE_TFs (decreases TEPIC TF affinities for TFs not found in TGENE at target gene)".
    public String tgen_consensus_calc="INCREASE_TGENE_TFS";
    //#[OPT]: if {OPT-SELF-REG} TGENE consensus is used please specify the writing of the Mitochondrial DNA chromosome in Peak Files, default: MT
    public String tgen_mt_writing="MT";

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

    /*######################
    ####PLOT parameters###
    ######################*/
    //#[OPT]: thresholds for coefficent plot and overall plots, for each threshold it creates a plot for each timepoint and an overall plot for each HistoneMod
    //#e.g. 0.1 means it creates a plot of the coefficient range ([-1;1]), it uses all TFs of [-1,-0.1] and [0.1,1]
    //#default: 0.1;0.2;0.3;0.4;0.5;0.6
    public List<Double> plot_th_coefficient= new ArrayList<Double>(Arrays.asList(new Double[]{0.1, 0.2, 0.3, 0.4, 0.5, 0.6}));
    //#[OPT]: in how many TPs should a TF be found, default: 1
    public int plot_cutoff_tps=1;
    //#[OPT]: in how many HMs should a TF be found, default: 1
    public int plot_cutoff_hms=1;
    //#[OPT]: minimum gene counts, default: 100
    public int plot_cutoff_gcs=100;
    //#[OPT]: minimum TPM, default: 0.0
    public double plot_cutoff_tpms=0.0;
    //#[OPT]: top k target genes for TFs, default: 30
    public int plot_top_k_genes=30;
    //#[REQ]: mann-whitneyU pvalue cutoff (usually 0.05 or 0.01), default 0.01
    public double plot_mann_whitneyU_pvalue_cutoff=0.01;

    /*#########################
      ####WEBSITE parameters###
      #########################*/
    //#[OPT]: list of TFs which you are interested in - is only used to search fast for known TFs in the results, it does not affect results
    //#e.g. website_interesting_tfs="STAT3;GATA3;NFIB"
    public HashSet<String> website_interesting_tfs=new HashSet<>();

    /*############################
      ####ChIP-atlas parameters###
      ############################*/
    public boolean chip_atlas_activated_chip_atlas =false;
    //before the restructuring of CHIP-Atlas
    //public String chip_atlas_url_to_list="https://dbarchive.biosciencedbc.jp/data/chip-atlas/LATEST/chip_atlas_file_list.zip";
    public String chip_atlas_url_to_list="http://togodb.biosciencedbc.jp/togodb/release/chip_atlas_file_list.csv";
    public String chip_atlas_column_gene_version=".*genome.*assembly.*";
    public String chip_atlas_column_antigen_class=".*antigen.*class.*";
    public String chip_atlas_column_antigen=".*antigen.*";
    public String chip_atlas_column_cell_type_class=".*cell.*type.*class.*";
    public String chip_atlas_column_url=".*file.*url.*";
    //##if no validation TF ChIP-seq data is available, you can provide the chip atlas genome version and the tissue type, to automatically get all available ChIP atlas data for the TFs
    //#[OPT]: chip atlas genome version (e.g. hg19, hg38, mm9, mm10, ...) look at: https://chip-atlas.org/peak_browser
    public String chip_atlas_genome_version="";
    //#[OPT]: chip atlas tissue type (e.g. lung, ...), look at: https://chip-atlas.org/peak_browser
    public String chip_atlas_tissue_type="";

    /*########################
        ####IGV parameters###
        #########################*/
    //#[OPT]: path to igv.sh, if not set no images of the top target genes of the discounted cumulative gain TFs will be taken
    public String igv_path_to_igv="";
    //#[OPT]: path to tf chip_seq data, of which you want to take images for target genes to see peaks
    public String igv_path_to_tf_chip_seq="";
    //#[OPT] port number of igv, default = 60151
    public int igv_port_number=60151;
    //#[OPT] igv species name for reference genome
    public String igv_species_ref_genome="";


    /*##################################
    ####after run analysis parameters###
    ####################################*/
    //these parameters are only used in ~/COM2POSE/src/further_analysis_tools/TPM_GC_Filter_Analysis
    //[REQ]: root-run-directories where all runs are in different folder including a working_dir folder in each folder
    public String tpm_gc_filter_analysis_working_dir = "";
    //[REQ]: TF-list file directory
    public String tpm_gc_filter_analysis_tf_list = "";
    //[OPT]: count zeros in means
    public boolean tpm_gc_filter_analysis_count_zeros=true;

    //dir_names
    public String tpm_gc_filter_analysis_folder_directory_name = "A1_TPM_GC_FILTER_ANALYSIS";
    public String tpm_gc_filter_analysis_folder_directory_name_data = "01_DATA";
    public String tpm_gc_filter_analysis_folder_directory_name_RScripts = "02_RScripts";
    public String tpm_gc_filter_analysis_folder_directory_name_plots = "03_PLOTS";

    public String tpm_gc_filter_analysis_suffix_data = "data.txt";

    /*##################################
    ######conversion helpers############
    ####################################*/

    public String rna_seq_conversion_input="";
    public String rna_seq_conversion_biomart_species_name="";
    public String rna_seq_conversion_biomart_column="";



}
