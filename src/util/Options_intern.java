package util;

public class Options_intern
{
    /**
     * COM2POSE OPTIONS intern
     */
    public String config_data_path="";

    /*
    ########################
    ##COM2POSE parameters###
    ########################*/
    //#[REQ]: working directory where COM2POSE can create, remove and edit files
    public String com2pose_working_directory="";
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
    //#[REQ]: path to position specific energymatrix used for TRAP (different matrices can be found in ~/COM2POSE/ext/TEPIC/TEPIC/PWMs/2.1)
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
    public int tepic_tpm_cutoff=-1;
    //#[OPT]: path to input file ensg to gene symbol file, required for TPM filter
    public String tepic_ensg_symbol="";
   /*######################
     ##DYNAMITE parameters#
     ######################*/



    /**
     * TEPIC OPTIONS

    public String input_dir_peaks = "";
    public String input_reference_genome = "";
    public String output_dir = "";
    public String gene_annotation = "";

    /**
     * OPTIONAL OPTIONS TEPIC

    public int window = 50000;
    public int cores = 1;

    /**
     * FLAGS TEPIC

    public boolean run_all_subdirectories = false;
     */
}
