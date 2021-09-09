package tcga_atac_seq_helpers;

import java.io.*;
import java.util.*;

public class rna_seq_unpack_sort {

    public static String input_directory_biogrid_entrez= "/home/markus/Dropbox/UNI/Promotion_Projekte/PhD_work/H_ATAC_SEQ/PCA/olga_analysis/biogrid.human.entrez.tsv";

    public static String action ="";
    public static String only_execute_samples ="";
    public static String execute_option ="";

    public static void main(String[] args) throws Exception {

        //EXECUTE = execute all sorting, formatting and pca
        //REMOVE = removes all files created by this pipeline for a new try afterwards
        //ANALYSE = gives back all samples_which are not 01A (primary tumor RNA-seq)
        action ="EXECUTE";


        //; seperated list of which sets should be executed
        //if empty all will be executed
        only_execute_samples ="BRCA";

        //empty all will be executed
        //01A excludes all other than 01A
        //01B excludes all other than 01B
        //and so on
        execute_option="01A";



        /**
         * TGCT
         */
        atac_seq_download_cancer_blocks block_tgct = new atac_seq_download_cancer_blocks();
        {
            block_tgct.name="RNA_SEQ: TGCT";
            block_tgct.input_dir_root="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/TGCT/RNA_SEQ/01_raw/gdc_download_20210810_132716.957271";
            block_tgct.input_clinical="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/TGCT/METADATA/RNA_SEQ/clinical.cart.2021-08-10/clinical.tsv";
            block_tgct.input_gdc_sample_sheet = "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/TGCT/METADATA/RNA_SEQ/gdc_sample_sheet.2021-08-10.tsv";
            block_tgct.output_directory= "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/TGCT";

            //groups to which to sort to
            block_tgct.groups = new HashSet<>();
            block_tgct.groups.add("Seminoma");
            //name for other samples, which do not fit to groups
            block_tgct.other_samples = "NonSeminoma";
            block_tgct.number_groups =2;
            execute_all(block_tgct);
        }


        /**
         * BLCA
         */
        atac_seq_download_cancer_blocks block_BLCA = new atac_seq_download_cancer_blocks();
        {
            block_BLCA.name="RNA_SEQ: BLCA";
            block_BLCA.input_dir_root="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/03_BLCA/RNA_SEQ/gdc_download_20210811_064600.308234";
            block_BLCA.input_clinical="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/03_BLCA/META_DATA/RNA_SEQ/clinical.cart.2021-08-11/clinical.tsv";
            block_BLCA.input_gdc_sample_sheet = "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/03_BLCA/META_DATA/RNA_SEQ/gdc_sample_sheet.2021-08-11.tsv";
            block_BLCA.output_directory= "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/03_BLCA";

            //groups to which to sort to
            block_BLCA.groups = new HashSet<>();
            block_BLCA.groups.add("Papillary");
            //name for other samples, which do not fit to groups
            block_BLCA.other_samples = "NonPapillary";
            block_BLCA.number_groups=2;

            //EXECUTE ALL HAS TO BE UNCOMMENTED!
            execute_all(block_BLCA);
        }

        /**
         * BRCA
         */
        atac_seq_download_cancer_blocks block_BRCA = new atac_seq_download_cancer_blocks();
        {
            block_BRCA.name="RNA_SEQ: BRCA";
            block_BRCA.input_dir_root="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/03_BRCA/RNA_SEQ/01_raw/gdc_download_20210806_054044.782364";
            block_BRCA.input_clinical="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/03_BRCA/META_DATA/RNA_SEQ/clinical.cart.2021-08-06/clinical.tsv";
            block_BRCA.input_gdc_sample_sheet = "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/03_BRCA/META_DATA/RNA_SEQ/gdc_sample_sheet.2021-08-06.tsv";
            block_BRCA.output_directory= "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/03_BRCA";

            //groups to which to sort to
            block_BRCA.groups = new HashSet<>();
            block_BRCA.groups.add("Lobular");
            block_BRCA.groups.add("Infiltrating");
            block_BRCA.groups.add("Metaplastic");
            //name for other samples, which do not fit to groups
            block_BRCA.other_samples = "Others";
            block_BRCA.number_groups=4;
            execute_all(block_BRCA);
        }

        /**
         * COAD
         */
        atac_seq_download_cancer_blocks block_COAD = new atac_seq_download_cancer_blocks();
        {
            block_COAD.name="RNA_SEQ: COAD";
            block_COAD.input_dir_root="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/04_COAD/RNA_SEQ/01_raw/gdc_download_20210811_064850.270781";
            block_COAD.input_clinical="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/04_COAD/METADATA/RNA_SEQ/clinical.cart.2021-08-11/clinical.tsv";
            block_COAD.input_gdc_sample_sheet = "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/04_COAD/METADATA/RNA_SEQ/gdc_sample_sheet.2021-08-11.tsv";
            block_COAD.output_directory= "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/04_COAD";

            //groups to which to sort to
            block_COAD.groups = new HashSet<>();
            block_COAD.groups.add("Adenocarcinoma");
            block_COAD.other_samples = "NonAdenocarcinoma";
            block_COAD.number_groups=2;

            execute_all(block_COAD);
        }

        /**
         * KIRP
         */
        atac_seq_download_cancer_blocks block_KIRP = new atac_seq_download_cancer_blocks();
        {
            block_KIRP.name="RNA_SEQ: KIRP";
            block_KIRP.input_dir_root="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/05_KIRP/RNA_SEQ/01_raw/gdc_download_20210811_065137.607520";
            block_KIRP.input_clinical="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/05_KIRP/METADATA/RNA_SEQ/clinical.cart.2021-08-11/clinical.tsv";
            block_KIRP.input_gdc_sample_sheet = "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/05_KIRP/METADATA/RNA_SEQ/gdc_sample_sheet.2021-08-11.tsv";
            block_KIRP.output_directory= "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/05_KIRP";

            //groups to which to sort to
            block_KIRP.groups = new HashSet<>();
            block_KIRP.groups.add("Papillary");
            block_KIRP.other_samples = "NonPapillary";
            block_KIRP.number_groups=-1;

            execute_all(block_KIRP);
        }

        /**
         * PRAD
         */
        atac_seq_download_cancer_blocks block_PRAD = new atac_seq_download_cancer_blocks();
        {
            block_PRAD.name="RNA_SEQ: PRAD";
            block_PRAD.input_dir_root="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/06_PRAD/RNA_SEQ/01_raw/gdc_download_20210811_065355.069019";
            block_PRAD.input_clinical="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/06_PRAD/METADATA/RNA_SEQ/clinical.cart.2021-08-11/clinical.tsv";
            block_PRAD.input_gdc_sample_sheet = "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/06_PRAD/METADATA/RNA_SEQ/gdc_sample_sheet.2021-08-11.tsv";
            block_PRAD.output_directory= "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/06_PRAD";

            //groups to which to sort to
            block_PRAD.groups = new HashSet<>();
            block_PRAD.groups.add("Adenocarcinoma");
            block_PRAD.other_samples = "Infiltrating";
            block_PRAD.number_groups=-1;

            execute_all(block_PRAD);
        }

        /**
         * LUAD
         */
        atac_seq_download_cancer_blocks block_LUAD = new atac_seq_download_cancer_blocks();
        {
            block_LUAD.name="RNA_SEQ: LUAD";
            block_LUAD.input_dir_root="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/07_LUAD/RNA_SEQ/01_raw/gdc_download_20210811_071055.539794";
            block_LUAD.input_clinical="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/07_LUAD/METADATA/RNA_SEQ/clinical.cart.2021-08-11/clinical.tsv";
            block_LUAD.input_gdc_sample_sheet = "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/07_LUAD/METADATA/RNA_SEQ/gdc_sample_sheet.2021-08-11.tsv";
            block_LUAD.output_directory= "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/07_LUAD";

            //groups to which to sort to
            block_LUAD.groups = new HashSet<>();
            block_LUAD.groups.add("Adenocarcinoma");
            block_LUAD.groups.add("Acinar");
            block_LUAD.groups.add("Mucinous");
            block_LUAD.other_samples = "PapillaryAdenocarcinoma";
            block_LUAD.number_groups=-1;

            execute_all(block_LUAD);
        }

        /**
         * STAD
         */
        atac_seq_download_cancer_blocks block_STAD = new atac_seq_download_cancer_blocks();
        {
            block_STAD.name="RNA_SEQ: STAD";
            block_STAD.input_dir_root="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/08_STAD/RNA_SEQ/01_raw/gdc_download_20210811_065952.752534";
            block_STAD.input_clinical="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/08_STAD/METADATA/RNA_SEQ/clinical.cart.2021-08-11/clinical.tsv";
            block_STAD.input_gdc_sample_sheet = "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/08_STAD/METADATA/RNA_SEQ/gdc_sample_sheet.2021-08-11.tsv";
            block_STAD.output_directory= "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/08_STAD";

            //groups to which to sort to
            block_STAD.groups = new HashSet<>();
            block_STAD.groups.add("Adenocarcinoma");
            block_STAD.groups.add("Mucinous");
            block_STAD.groups.add("Tubular");
            block_STAD.other_samples = "NonAdenocarcinoma";
            block_STAD.number_groups=-1;

            execute_all(block_STAD);
        }

        /**
         * ESCA
         */
        atac_seq_download_cancer_blocks block_ESCA = new atac_seq_download_cancer_blocks();
        {
            block_ESCA.name="RNA_SEQ: ESCA";
            block_ESCA.input_dir_root="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/09_ESCA/RNA_SEQ/01_raw/gdc_download_20210811_070227.439914";
            block_ESCA.input_clinical="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/09_ESCA/METADATA/RNA_SEQ/clinical.cart.2021-08-11/clinical.tsv";
            block_ESCA.input_gdc_sample_sheet = "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/09_ESCA/METADATA/RNA_SEQ/gdc_sample_sheet.2021-08-11.tsv";
            block_ESCA.output_directory= "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/09_ESCA";

            //groups to which to sort to
            block_ESCA.groups = new HashSet<>();
            block_ESCA.groups.add("Squamous");
            block_ESCA.other_samples = "Adenocarcinoma";
            block_ESCA.number_groups=-1;

            execute_all(block_ESCA);
        }

        /**
         * LIHC
         */
        atac_seq_download_cancer_blocks block_LIHC = new atac_seq_download_cancer_blocks();
        {
            block_LIHC.name="RNA_SEQ: LIHC";
            block_LIHC.input_dir_root="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/10_LIHC/RNA_SEQ/01_raw/gdc_download_20210811_071408.358464";
            block_LIHC.input_clinical="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/10_LIHC/METADATA/RNA_SEQ/clinical.cart.2021-08-11/clinical.tsv";
            block_LIHC.input_gdc_sample_sheet = "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/10_LIHC/METADATA/RNA_SEQ/gdc_sample_sheet.2021-08-11.tsv";
            block_LIHC.output_directory= "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/10_LIHC";

            //groups to which to sort to
            block_LIHC.groups = new HashSet<>();
            block_LIHC.groups.add("Hepatocellular");
            block_LIHC.other_samples = "NonHepatocellular";
            block_LIHC.number_groups=-1;

            execute_all(block_LIHC);
        }

        /**
         * KIRC
         */
        atac_seq_download_cancer_blocks block_KIRC = new atac_seq_download_cancer_blocks();
        {
            block_KIRC.name="RNA_SEQ: KIRC";
            block_KIRC.input_dir_root="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/11_KIRC/RNA_SEQ/01_raw/gdc_download_20210811_071630.749633";
            block_KIRC.input_clinical="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/11_KIRC/METADATA/RNA_SEQ/clinical.cart.2021-08-11/clinical.tsv";
            block_KIRC.input_gdc_sample_sheet = "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/11_KIRC/METADATA/RNA_SEQ/gdc_sample_sheet.2021-08-11.tsv";
            block_KIRC.output_directory= "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/11_KIRC";

            //groups to which to sort to
            block_KIRC.groups = new HashSet<>();
            block_KIRC.groups.add("Clear");
            block_KIRC.other_samples = "Renal,";
            block_KIRC.number_groups=-1;

            execute_all(block_KIRC);
        }

        /**
         * LUSC
         */
        atac_seq_download_cancer_blocks block_LUSC = new atac_seq_download_cancer_blocks();
        {
            block_LUSC.name="RNA_SEQ: LUSC";
            block_LUSC.input_dir_root="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/12_LUSC/RNA_SEQ/01_raw/gdc_download_20210811_071834.296668";
            block_LUSC.input_clinical="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/12_LUSC/METADATA/RNA_SEQ/clinical.cart.2021-08-11/clinical.tsv";
            block_LUSC.input_gdc_sample_sheet = "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/12_LUSC/METADATA/RNA_SEQ/gdc_sample_sheet.2021-08-11.tsv";
            block_LUSC.output_directory= "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/12_LUSC";

            //groups to which to sort to
            block_LUSC.groups = new HashSet<>();
            block_LUSC.groups.add("keratinizing");
            block_LUSC.other_samples = "Squamous";
            block_LUSC.number_groups=-1;

            execute_all(block_LUSC);
        }

        /**
         * THCA
         */
        atac_seq_download_cancer_blocks block_THCA = new atac_seq_download_cancer_blocks();
        {
            block_THCA.name="RNA_SEQ: THCA";
            block_THCA.input_dir_root="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/13_THCA/RNA_SEQ/01_raw/gdc_download_20210811_072143.524926";
            block_THCA.input_clinical="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/13_THCA/METADATA/RNA_SEQ/clinical.cart.2021-08-11/clinical.tsv";
            block_THCA.input_gdc_sample_sheet = "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/13_THCA/METADATA/RNA_SEQ/gdc_sample_sheet.2021-08-11.tsv";
            block_THCA.output_directory= "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/13_THCA";

            //groups to which to sort to
            block_THCA.groups = new HashSet<>();
            block_THCA.groups.add("Papillary");
            block_THCA.groups.add("Oxyphilic");
            block_THCA.other_samples = "Follicular";
            block_THCA.number_groups=-1;

            execute_all(block_THCA);
        }

        /**
         * LGG
         */
        atac_seq_download_cancer_blocks block_LGG = new atac_seq_download_cancer_blocks();
        {
            block_LGG.name="RNA_SEQ: LGG";
            block_LGG.input_dir_root="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/14_LGG/RNA_SEQ/01_raw/gdc_download_20210811_072408.947445";
            block_LGG.input_clinical="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/14_LGG/METADATA/RNA_SEQ/clinical.cart.2021-08-11/clinical.tsv";
            block_LGG.input_gdc_sample_sheet = "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/14_LGG/METADATA/RNA_SEQ/gdc_sample_sheet.2021-08-11.tsv";
            block_LGG.output_directory= "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/14_LGG";

            //groups to which to sort to
            block_LGG.groups = new HashSet<>();
            block_LGG.groups.add("Oligodendroglioma");
            block_LGG.groups.add("Astrocytoma");
            block_LGG.other_samples = "Mixed";
            block_LGG.number_groups=-1;

            execute_all(block_LGG);
        }

        /**
         * SKCM
         */
        atac_seq_download_cancer_blocks block_SKCM = new atac_seq_download_cancer_blocks();
        {
            block_SKCM.name="RNA_SEQ: SKCM";
            block_SKCM.input_dir_root="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/15_SKCM/RNA_SEQ/01_raw/gdc_download_20210811_072645.075074";
            block_SKCM.input_clinical="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/15_SKCM/METADATA/RNA_SEQ/clinical.cart.2021-08-11/clinical.tsv";
            block_SKCM.input_gdc_sample_sheet = "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/15_SKCM/METADATA/RNA_SEQ/gdc_sample_sheet.2021-08-11.tsv";
            block_SKCM.output_directory= "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/15_SKCM";

            //groups to which to sort to
            block_SKCM.groups = new HashSet<>();
            block_SKCM.groups.add("Malignant");
            block_SKCM.other_samples = "Epithelioid";
            block_SKCM.number_groups=-1;

            execute_all(block_SKCM);
        }

        /**
         * UCEC
         */
        atac_seq_download_cancer_blocks block_UCEC = new atac_seq_download_cancer_blocks();
        {
            block_UCEC.name="RNA_SEQ: UCEC";
            block_UCEC.input_dir_root="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/16_UCEC/RNA_SEQ/01_raw/gdc_download_20210811_072934.082004";
            block_UCEC.input_clinical="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/16_UCEC/METADATA/RNA_SEQ/clinical.cart.2021-08-11/clinical.tsv";
            block_UCEC.input_gdc_sample_sheet = "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/16_UCEC/METADATA/RNA_SEQ/gdc_sample_sheet.2021-08-11.tsv";
            block_UCEC.output_directory= "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/16_UCEC";

            //groups to which to sort to
            block_UCEC.groups = new HashSet<>();
            block_UCEC.groups.add("Endometrioid");
            block_UCEC.other_samples = "SerousCystadenocarcinoma";
            block_UCEC.number_groups=-1;

            execute_all(block_UCEC);
        }

        /**
         * ACC
         */
        atac_seq_download_cancer_blocks block_ACC = new atac_seq_download_cancer_blocks();
        {
            block_ACC.name="RNA_SEQ: ACC";
            block_ACC.input_dir_root="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/17_ACC/RNA_SEQ/01_raw/gdc_download_20210811_073135.623128";
            block_ACC.input_clinical="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/17_ACC/METADATA/RNA_SEQ/clinical.cart.2021-08-11/clinical.tsv";
            block_ACC.input_gdc_sample_sheet = "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/17_ACC/METADATA/RNA_SEQ/gdc_sample_sheet.2021-08-11.tsv";
            block_ACC.output_directory= "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/17_ACC";

            //groups to which to sort to
            block_ACC.groups = new HashSet<>();
            block_ACC.groups.add("Adrenal");
            block_ACC.other_samples = "NonAdrenal";
            block_ACC.number_groups=-1;

            execute_all(block_ACC);
        }

        /**
         * GBM IS MISSING RNA-SEQ DATA
         */

        /**
         * HNSC
         */
        atac_seq_download_cancer_blocks block_HNSC = new atac_seq_download_cancer_blocks();
        {
            block_HNSC.name="RNA_SEQ: HNSC";
            block_HNSC.input_dir_root="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/19_HNSC/RNA_SEQ/01_raw/gdc_download_20210811_073652.791948";
            block_HNSC.input_clinical="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/19_HNSC/METADATA/RNA_SEQ/clinical.cart.2021-08-11/clinical.tsv";
            block_HNSC.input_gdc_sample_sheet = "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/19_HNSC/METADATA/RNA_SEQ/gdc_sample_sheet.2021-08-11.tsv";
            block_HNSC.output_directory= "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/19_HNSC";

            //groups to which to sort to
            block_HNSC.groups = new HashSet<>();
            block_HNSC.groups.add("Basaloid");
            block_HNSC.groups.add("keratinizing");
            block_HNSC.other_samples = "Squamous";
            block_HNSC.number_groups=-1;

            execute_all(block_HNSC);
        }

        /**
         * PCPG
         */
        atac_seq_download_cancer_blocks block_PCPG = new atac_seq_download_cancer_blocks();
        {
            block_PCPG.name="RNA_SEQ: PCPG";
            block_PCPG.input_dir_root="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/20_PCPG/RNA_SEQ/01_raw/gdc_download_20210811_073906.011732";
            block_PCPG.input_clinical="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/20_PCPG/METADATA/RNA_SEQ/clinical.cart.2021-08-11/clinical.tsv";
            block_PCPG.input_gdc_sample_sheet = "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/20_PCPG/METADATA/RNA_SEQ/gdc_sample_sheet.2021-08-11.tsv";
            block_PCPG.output_directory= "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/20_PCPG";

            //groups to which to sort to
            block_PCPG.groups = new HashSet<>();
            block_PCPG.groups.add("malignant");
            block_PCPG.other_samples = "NonMalignant";
            block_PCPG.number_groups=-1;

            execute_all(block_PCPG);
        }

        /**
         * MESO
         */
        atac_seq_download_cancer_blocks block_MESO = new atac_seq_download_cancer_blocks();
        {
            block_MESO.name="RNA_SEQ: MESO";
            block_MESO.input_dir_root="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/21_MESO/RNA_SEQ/01_raw/gdc_download_20210811_074127.501831";
            block_MESO.input_clinical="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/21_MESO/METADATA/RNA_SEQ/clinical.cart.2021-08-11/clinical.tsv";
            block_MESO.input_gdc_sample_sheet = "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/21_MESO/METADATA/RNA_SEQ/gdc_sample_sheet.2021-08-11.tsv";
            block_MESO.output_directory= "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/21_MESO";

            //groups to which to sort to
            block_MESO.groups = new HashSet<>();
            block_MESO.groups.add("biphasic");
            block_MESO.groups.add("Epithelioid");
            block_MESO.other_samples = "Mesothelioma";
            block_MESO.number_groups=-1;

            execute_all(block_MESO);
        }

        /**
         * CHOL
         */
        atac_seq_download_cancer_blocks block_CHOL = new atac_seq_download_cancer_blocks();
        {
            block_CHOL.name="RNA_SEQ: CHOL";
            block_CHOL.input_dir_root="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/22_CHOL/RNA_SEQ/01_raw/gdc_download_20210811_074351.486592";
            block_CHOL.input_clinical="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/22_CHOL/METADATA/RNA_SEQ/clinical.cart.2021-08-11/clinical.tsv";
            block_CHOL.input_gdc_sample_sheet = "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/22_CHOL/METADATA/RNA_SEQ/gdc_sample_sheet.2021-08-11.tsv";
            block_CHOL.output_directory= "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/22_CHOL";

            //groups to which to sort to
            block_CHOL.groups = new HashSet<>();
            block_CHOL.groups.add("Cholangiocarcinoma");
            block_CHOL.other_samples = "NonCholangiocarcinoma";
            block_CHOL.number_groups=-1;

            execute_all(block_CHOL);
        }

        /**
         * CESC
         */
        atac_seq_download_cancer_blocks block_CESC = new atac_seq_download_cancer_blocks();
        {
            block_CESC.name="RNA_SEQ: CESC";
            block_CESC.input_dir_root="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/23_CESC/RNA_SEQ/01_raw/gdc_download_20210811_074548.241451";
            block_CESC.input_clinical="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/23_CESC/METADATA/RNA_SEQ/clinical.cart.2021-08-11/clinical.tsv";
            block_CESC.input_gdc_sample_sheet = "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/23_CESC/METADATA/RNA_SEQ/gdc_sample_sheet.2021-08-11.tsv";
            block_CESC.output_directory= "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/23_CESC";

            //groups to which to sort to
            block_CESC.groups = new HashSet<>();
            block_CESC.groups.add("Squamous");
            block_CESC.other_samples = "NonSquamous";
            block_CESC.number_groups=-1;

            execute_all(block_CESC);
        }

        atac_seq_download_cancer_blocks block_EXAMPLE = new atac_seq_download_cancer_blocks();
        {
            block_EXAMPLE.name="RNA_SEQ: XXXX";
            block_EXAMPLE.input_dir_root="";
            block_EXAMPLE.input_clinical="/clinical.tsv";
            block_EXAMPLE.input_gdc_sample_sheet = "";
            block_EXAMPLE.output_directory= "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/";

            //groups to which to sort to
            block_EXAMPLE.groups = new HashSet<>();
            block_EXAMPLE.groups.add("XXX");
            //name for other samples, which do not fit to groups
            block_EXAMPLE.other_samples = "Others";
            block_EXAMPLE.number_groups=-1;

            //execute_all(block_EXAMPLE);
        }




    }

    private static void execute_all(atac_seq_download_cancer_blocks blocks) throws Exception {

        if(!only_execute_samples.equals(""))
        {
           String[] split_only_execute =  only_execute_samples.split(";");
           Set<String> execute_these = new HashSet<String>(Arrays.asList(split_only_execute));


           String[] split_blocks_name=blocks.name.split(":");
           String current_name= split_blocks_name[1].trim();

           if(!execute_these.contains(current_name))
           {
               return;
           }
        }

        blocks.number_groups=blocks.groups.size()+1;

        String input_sample_sheet=blocks.input_gdc_sample_sheet;
        String input_clinical = blocks.input_clinical;
        String output_sorted=blocks.output_directory;
        String input_root=blocks.input_dir_root;

        File f_output_sorted_root = new File(output_sorted);
        File f_output_rnaseq = new File(f_output_sorted_root.getAbsolutePath()+File.separator+"RNA_SEQ");
        f_output_rnaseq.mkdirs();
        File f_output_rnaseq_sorted = new File(f_output_rnaseq.getAbsolutePath()+File.separator+"02_ordered");
        f_output_rnaseq_sorted.mkdirs();
        File f_output_rnaseq_formatted = new File(f_output_rnaseq.getAbsolutePath()+File.separator+"03_formatted");
        f_output_rnaseq_formatted.mkdirs();

        File f_output_pca = new File(f_output_rnaseq.getAbsolutePath()+File.separator+"04_PCA_UMAP");
        f_output_pca.mkdir();

        File f_output_mapping = new File(f_output_rnaseq.getAbsolutePath()+File.separator+"mapping_samples.csv");

        if(action.equals("ANALYSE"))
        {
            if(!f_output_mapping.exists())
            {
                System.out.println("[ERROR] mapping does not exist - action EXECUTE has to be done first!" + blocks.name);
                return;
            }

            HashMap<String,HashSet<String>> ident_name = analyse_it(f_output_mapping,"SAMPLE_NAME");


            System.out.println("X");



        }else if(action.equals("REMOVE"))
        {

            String command_remove_ordered = "rm -r " + f_output_rnaseq_sorted.getAbsolutePath();
            String command_remove_formatted = "rm -r " + f_output_rnaseq_formatted.getAbsolutePath()+File.separator;
            String command_remove_pca = "rm -r " + f_output_pca.getAbsolutePath();

            Process child_ordered = Runtime.getRuntime().exec(command_remove_ordered);
            int code_ordered = child_ordered.waitFor();
            switch (code_ordered){
                case 0:
                    break;
                case 1:
                    String message = child_ordered.getErrorStream().toString();
                    //throw new Exception(message);
            }

            Process child_formatted = Runtime.getRuntime().exec(command_remove_formatted);
            int code_formatted = child_formatted.waitFor();
            switch (code_formatted){
                case 0:
                    break;
                case 1:
                    String message = child_formatted.getErrorStream().toString();
                    //throw new Exception(message);
            }

            Process child_remove = Runtime.getRuntime().exec(command_remove_pca);
            int code_remove = child_remove.waitFor();
            switch (code_remove){
                case 0:
                    break;
                case 1:
                    String message = child_remove.getErrorStream().toString();
                    //throw new Exception(message);
            }

            if(code_formatted==0 && code_remove==0 && code_ordered==0)
            {
                System.out.println("[SUCCESS] removed all from "+ blocks.name);
            }
            else
            {
                System.out.println("[ERROR] could not remove all from "+ blocks.name);
            }


            return;
        }else if(action.equals("EXECUTE"))
        {

            File f_input_root = new File(input_root);

            System.out.println("PERFORM ON"+ blocks.name);

            HashMap<String,HashSet<String>> ident_folder = analyse_it(f_output_mapping, "FOLDER_NAME");

            boolean unpack_successfull = unpack(f_input_root);
            if(!unpack_successfull)
            {
                System.out.println("[ERROR] " + blocks.name + " did not unpack successfull... aborting ...");
                return;
            }

            File f_input_clinical = new File(input_clinical);
            File f_input_sample_sheet = new File(input_sample_sheet);
            if(!f_input_sample_sheet.exists())
            {
                System.out.println("X");
            }

            sort(f_input_root,f_output_sorted_root,f_output_rnaseq_sorted,blocks.groups,blocks.other_samples, f_input_clinical, f_input_sample_sheet, f_output_mapping, ident_folder, f_output_rnaseq_sorted, f_output_rnaseq_formatted, f_output_pca);

            //format_uneven_ensgs(f_output_rnaseq_sorted,f_output_rnaseq_formatted);

            format_even_ensgs(f_output_rnaseq_sorted,f_output_rnaseq_formatted);

            //perform PCA
            boolean pca_successfull = perform_pca(f_output_rnaseq_sorted, f_output_pca, blocks);

            if(!pca_successfull)
            {
                System.out.println("[ERROR] " + blocks.name + " pca error. please try manually");
                return;
            }

        }
        else
        {
            System.out.println("NO VALID ACTION CHOSEN: EXECUTE OR REMOVE");
        }




    }

    /**
     *
     * @param f_output_mapping mapping file
     * @param which_ident SAMPLE_NAME = my own created sample name; FOLDER_NAME of raw data; FILE_NAME = htseq count file name
     * @return sample_identifiers to sample names
     * @throws IOException
     */
    private static HashMap<String,HashSet<String>> analyse_it(File f_output_mapping, String which_ident) throws IOException {

        HashMap<String,HashSet<String>> ident_name = new HashMap<>();

        if(!f_output_mapping.exists())
        {
            return ident_name;
        }

        BufferedReader br_mapping = new BufferedReader(new FileReader(f_output_mapping));
        String line = br_mapping.readLine();
        String[] split_header = line.split("\t");
        while((line=br_mapping.readLine())!=null)
        {

            String[] split = line.split("\t");

            String sample_name="";
            if(which_ident.equals("SAMPLE_NAME"))
            {
                sample_name = split[0];
            }
            else if(which_ident.equals("FOLDER_NAME"))
            {
                sample_name = split[2];
            }
            else if(which_ident.equals("FILE_NAME"))
            {
                sample_name = split[4];
            }

            String identifier = split[7].split("-")[3];

            if(!ident_name.containsKey(identifier))
            {
                HashSet<String> x = new HashSet<>();
                x.add(sample_name);

                ident_name.put(identifier,x);
            }
            else
            {
                ident_name.get(identifier).add(sample_name);
            }
        }
        br_mapping.close();

        return ident_name;

    }

    private static boolean perform_pca(File f_output_rnaseq_sorted, File f_output_pca, atac_seq_download_cancer_blocks blocks) throws Exception {

        boolean correct = false;

        System.out.println("performing PCA analysis");

        File f_script_folder = new File(f_output_pca.getAbsolutePath()+File.separator+"01_script");
        f_script_folder.mkdir();
        File f_data_folder = new File(f_output_pca.getAbsolutePath()+File.separator+"02_pydata");
        f_data_folder.mkdir();

        File f_data_data = new File(f_data_folder.getAbsolutePath()+File.separator+"counts.csv");
        File f_data_meta = new File(f_data_folder.getAbsolutePath()+File.separator+"meta.csv");
        File f_data_out = new File(f_data_folder.getAbsolutePath()+File.separator+"out.csv");
        File f_data_norm = new File(f_data_folder.getAbsolutePath()+File.separator+"norm.csv");

        File f_script_file = new File(f_script_folder.getAbsolutePath()+File.separator+"pca_script.py");

        File f_plot_pca = new File(f_output_pca.getAbsolutePath()+File.separator+blocks.name+"_pca_plot.png");
        File f_plot_umap = new File(f_output_pca.getAbsolutePath()+File.separator+blocks.name+"_umap_plot.png");
        File f_plot_heatmap = new File(f_output_pca.getAbsolutePath()+File.separator+blocks.name+"_heatmap_plot.png");

        StringBuilder sb_python_script = new StringBuilder();

        sb_python_script.append("import pip\n" +
                "\n" +
                "def import_or_install(package):\n" +
                "    try:\n" +
                "        __import__(package)\n" +
                "    except ImportError:\n" +
                "        pip.main(['install', package])\n");

        sb_python_script.append("import_or_install(\"bicon\")\n" +
                "import_or_install(\"mygene\")\n" +
                "import_or_install(\"umap-learn\")\n" +
                "import_or_install(\"tqdm\")\n");

        sb_python_script.append("import pandas as pd\n" +
                "import matplotlib.pyplot as plt\n" +
                "\n" +
                "\n" +
                "path=\""+f_output_rnaseq_sorted.getAbsolutePath()+"\"\n" +
                "\n" +
                "import os\n" +
                "files = os.listdir(path)\n\n");

        sb_python_script.append("groups=[]\n" +
                "for file in files:\n" +
                "    if not file.endswith('.csv'):\n" +
                "        groups.append(file)\n" +
                "#print(groups)\n\n");
        sb_python_script.append("data=[]\n" +
                "meta=[]\n" +
                "for group in groups:\n" +
                "    group_path=path+\"/\"+group\n" +
                "    files = os.listdir(group_path)\n" +
                "    for file in files:\n" +
                "        file_path=group_path+\"/\"+file\n" +
                "        #print(\"file_path: \" + file_path)\n" +
                "        dat = pd.read_csv(file_path, sep = \"\\t\", header = None, index_col = 0)\n" +
                "        #print(dat)\n" +
                "        s_type=\"\".join(group)\n" +
                "        sample=group+file.split(\"_\")[2].split(\".\")[0]\n" +
                "        #print(\"s_type: \"+s_type)\n" +
                "        #print(\"sample: \" +sample)\n" +
                "        dat.columns = [sample]\n" +
                "        data.append(dat)\n" +
                "        meta.append((sample, s_type))\n\n");

        sb_python_script.append("data = pd.concat(data, axis = 1)\n" +
                "#data\n" +
                "meta = pd.DataFrame(meta,columns = [\"sample\", \"condition\"])\n" +
                "meta = meta.set_index(\"sample\")\n" +
                "data = data.T[data.index[:-5]].T\n" +
                "data.index.name = \"gene\"\n" +
                "data.to_csv(\""+f_data_data.getAbsolutePath()+"\")\n" +
                "meta.to_csv(\""+f_data_meta.getAbsolutePath()+"\")\n");

        sb_python_script.append("norm=(data-data.min())/(data.max()-data.min())\n");

        sb_python_script.append("from bicon import data_preprocessing\n" +
                "from bicon import BiCoN\n" +
                "from bicon import results_analysis\n" +
                "import mygene\n" +
                "mg = mygene.MyGeneInfo()\n\n");

        sb_python_script.append("index = norm.index\n" +
                "a_list = list(index)\n" +
                "#print(a_list)\n" +
                "\n" +
                "modified_list=[]\n" +
                "\n" +
                "for ensg in a_list:\n" +
                "    ensg=ensg.split(\".\")[0]\n" +
                "    modified_list.append(ensg)\n" +
                "\n" +
                "#print(modified_list)\n");

        sb_python_script.append("norm.index= modified_list\n\n");
        sb_python_script.append("out = mg.querymany(norm.index, scopes='ensemblgene', fields='entrezgene', species=\"human\", verbose=False, \n" +
                "                  entrezonly = True, as_dataframe = True)\n");
        sb_python_script.append("out.to_csv(\""+f_data_out.getAbsolutePath()+"\")\n");
        sb_python_script.append("new_genes = out[~out.entrezgene.isna()].index\n");
        sb_python_script.append("norm = norm.T[new_genes].T\n");
        sb_python_script.append("#out[~out.entrezgene.isna()]\n");
        sb_python_script.append("norm.index = out[~out.entrezgene.isna()].entrezgene\n");
        sb_python_script.append("norm.to_csv(\""+f_data_norm.getAbsolutePath()+"\")\n");
        sb_python_script.append("path_expr,path_net ='"+f_data_norm.getAbsolutePath()+"', '"+input_directory_biogrid_entrez+"'\n");
        sb_python_script.append("GE,G,labels, _= data_preprocessing(path_expr, path_net)\n");
        sb_python_script.append("L_g_min = 10\n" +
                "L_g_max = 15\n");
        sb_python_script.append("model = BiCoN(GE,G,L_g_min,L_g_max, )\n" +
                "solution,scores= model.run_search(show_plot = True)\n");
        sb_python_script.append("from sklearn.preprocessing import StandardScaler\n" +
                "from sklearn.decomposition import PCA\n");
        sb_python_script.append("import numpy as np\n" +
                "minimal = norm.min().min()\n" +
                "if minimal <= 0:\n" +
                "    norm += np.abs(minimal - 1)\n" +
                "norm = np.log2(norm)\n");
        sb_python_script.append("size = 2000\n" +
                "intersec_genes_list = norm.index\n" +
                "if size is not None:  # std selection\n" +
                "    std_genes = norm.std(axis=1)\n" +
                "    std_genes, intersec_genes_list = zip(*sorted(zip(std_genes, intersec_genes_list)))\n" +
                "    genes_for_expr = list(intersec_genes_list)[len(std_genes) - size:]\n" +
                "    intersec_genes = set(genes_for_expr)\n" +
                "\n" +
                "\n" +
                "norm = norm.loc[genes_for_expr]\n");
        sb_python_script.append("from scipy import stats\n" +
                "norm = pd.DataFrame((stats.zscore(norm.T)).T, columns=norm.columns, index=norm.index)\n");
        sb_python_script.append("pca = PCA(n_components="+blocks.number_groups+")\n" +
                "components = pca.fit_transform(norm.T)\n");

        sb_python_script.append("import seaborn as sns\n" +
                "fig, ax = plt.subplots(figsize=(10, 8))\n" +
                "\n" +
                "conditions = set(meta.condition)\n" +
                "colors = sns.color_palette(\"tab10\", len(conditions))\n" +
                "\n" +
                "c_lst = list(meta.condition)\n" +
                "for i,c in enumerate(conditions):\n" +
                "    idx = [i for i in range(len(c_lst)) if c_lst[i] == c]\n" +
                "    plt.scatter(components[idx,0],components[idx,1], color = colors[i], label = c)\n" +
                "plt.tight_layout()\n" +
                "plt.legend()\n" +
                "\n" +
                "#plt.show()\n");
        sb_python_script.append("plt.savefig(f'"+f_plot_pca.getAbsolutePath()+"')\n");

        sb_python_script.append("import umap.umap_ as umap\n");
        sb_python_script.append("reducer = umap.UMAP()\n" +
                "components = reducer.fit_transform(norm.T)\n");
        sb_python_script.append("fig, ax = plt.subplots(figsize=(10, 8))\n" +
                "\n" +
                "conditions = set(meta.condition)\n" +
                "colors = sns.color_palette(\"tab10\", len(conditions))\n" +
                "\n" +
                "c_lst = list(meta.condition)\n" +
                "for i,c in enumerate(conditions):\n" +
                "    idx = [i for i in range(len(c_lst)) if c_lst[i] == c]\n" +
                "    plt.scatter(components[idx,0],components[idx,1], color = colors[i], label = c)\n" +
                "plt.tight_layout()\n" +
                "plt.legend()\n");
        //save plot is missing
        sb_python_script.append("plt.savefig(f'"+f_plot_umap.getAbsolutePath()+"')\n");

        sb_python_script.append("grouping_p = meta\n" +
                "grouping_p.columns = [\"clusters\"]\n" +
                "species = grouping_p[\"clusters\"]\n" +
                "lut = {list(conditions)[i]:colors[i] for i in range(len(conditions))}\n" +
                "row_colors1 = species.map(lut)\n" +
                "\n" +
                "plt.rc('font', size=5)  # controls default text sizes\n" +
                "plt.rc('axes', titlesize=20)  # fontsize of the axes title\n" +
                "plt.rc('axes', labelsize=20)  # fontsize of the x and y labels\n" +
                "plt.rc('xtick', labelsize=20)  # fontsize of the tick labels\n" +
                "plt.rc('ytick', labelsize=10)  # fontsize of the tick labels\n" +
                "plt.rc('legend', fontsize=20)\n" +
                "\n" +
                "g = sns.clustermap(norm.T, row_colors=row_colors1, row_cluster=True, col_cluster=True,figsize=(15, 10),\n" +
                "                   cmap=\"Spectral\")\n" +
                "\n" +
                "\n" +
                "for key in lut:\n" +
                "    l = key\n" +
                "    c = lut[key]\n" +
                "    g.ax_col_dendrogram.bar(0, 0, color=c,\n" +
                "                            label=l, linewidth=0)\n" +
                "g.ax_col_dendrogram.legend(loc=\"upper center\", ncol=2, bbox_to_anchor=(0.72, 0.87),\n" +
                "                           borderaxespad=0.)\n" +
                "\n" +
                "ax = g.ax_heatmap\n" +
                "ax.set_xlabel(\"Genes\")\n" +
                "ax.set_ylabel(\"Patients\")\n" +
                "\n" +
                "#plt.show()\n");
        sb_python_script.append("plt.savefig(f'"+f_plot_heatmap.getAbsolutePath()+"')\n");


        BufferedWriter bw_script = new BufferedWriter(new FileWriter(f_script_file));
        bw_script.write(sb_python_script.toString());
        bw_script.close();

        //execute this script
        System.out.println("waiting...");

        String command_edited= "python3 " + f_script_file.getAbsolutePath();
        Process child = Runtime.getRuntime().exec(command_edited);
        int code = child.waitFor();
        switch (code){
            case 0:
                correct=true;
                break;
            case 1:
                correct=false;
                String message = child.getErrorStream().toString();
                //throw new Exception(message);
        }

        return correct;
    }

    private static void format_even_ensgs(File f_output_rnaseq_sorted, File f_output_rnaseq_formatted) throws IOException {

        HashMap<String,String> avail_ENSGs_files = new HashMap<>();

        for(File f_sorted : f_output_rnaseq_sorted.listFiles()) {
            if (f_sorted.isDirectory()) {
                for (File f_file : f_sorted.listFiles()) {
                    if (f_file.isFile()) {
                        BufferedReader br = new BufferedReader(new FileReader(f_file));
                        String line = "";
                        while((line=br.readLine())!=null)
                        {
                            if(line.startsWith("_"))
                                continue;
                            String[] split = line.split("\t");
                            avail_ENSGs_files.put(split[0].split("\\.")[0],split[0]);
                        }

                    }
                    break;
                }
            }
            break;
        }

        ArrayList<String> ensgs_ordered = new ArrayList<>(avail_ENSGs_files.keySet());
        Collections.sort(ensgs_ordered);
        HashMap<String,Integer> ensg_to_place = new HashMap<>();
        for(int i = 0; i < ensgs_ordered.size(); i++)
        {
            ensg_to_place.put(ensgs_ordered.get(i),i);
        }

        for(File f_sorted : f_output_rnaseq_sorted.listFiles())
        {
            if (f_sorted.isDirectory())
            {
                File f_output = new File(f_output_rnaseq_formatted.getAbsolutePath()+File.separator+f_sorted.getName());
                f_output.mkdirs();

                for (File f_file : f_sorted.listFiles())
                {
                    if (f_file.isFile())
                    {
                        HashMap<String,Boolean> already_got_ENSG = new HashMap<>();
                        HashMap<Integer,String> position_to_count = new HashMap<>();

                        for(String key_ensg : avail_ENSGs_files.keySet())
                        {
                            already_got_ENSG.put(key_ensg,false);
                        }

                        String[] rep_name_split= f_file.getName().split("_");
                        String replicate=rep_name_split[rep_name_split.length-1].split("\\.")[0];

                        String[] split_sorted = f_sorted.getName().split("_");
                        String groupname = "";
                        for(String k_split_sorted:split_sorted)
                        {
                            groupname+=k_split_sorted.substring(0, 1).toUpperCase() + k_split_sorted.substring(1);
                        }

                        String rep_name= groupname+"_Rep"+replicate+".txt";

                        File f_output_file = new File(f_output.getAbsolutePath()+File.separator+rep_name);

                        BufferedReader br = new BufferedReader(new FileReader(f_file));

                        String line = "";
                        while((line=br.readLine())!=null)
                        {
                            if(line.startsWith("_"))
                                continue;

                            String[] split = line.split("\t");
                            String ensg = split[0].split("\\.")[0];
                            if(already_got_ENSG.containsKey(ensg))
                            {
                                if(already_got_ENSG.get(ensg))
                                {
                                    System.out.println("DOUBLE " + ensg);
                                    continue;
                                }
                                else
                                {
                                    already_got_ENSG.put(ensg,true);

                                    int place_in_file = ensg_to_place.get(ensg);

                                    position_to_count.put(place_in_file,split[1]);
                                }
                            }
                        }
                        br.close();

                        BufferedWriter bw = new BufferedWriter(new FileWriter(f_output_file));
                        bw.write(rep_name);
                        bw.newLine();

                        for(int i = 0 ; i < ensgs_ordered.size(); i++)
                        {
                            if(!position_to_count.containsKey(i))
                            {
                                String ensg = ensgs_ordered.get(i);

                                System.out.println("not found: " + ensg);
                            }
                            bw.write(position_to_count.get(i));
                            bw.newLine();
                        }

                        bw.close();
                    }
                }
            }
        }

        File f_output_names_formatted = new File(f_output_rnaseq_formatted.getAbsolutePath()+File.separator+"names_formatted.txt");
        BufferedWriter bw_names = new BufferedWriter(new FileWriter(f_output_names_formatted));
        bw_names.write("Geneid");
        bw_names.newLine();
        for(int i = 0; i < ensgs_ordered.size();i++)
        {
            bw_names.write(ensgs_ordered.get(i));
            bw_names.newLine();
        }
        bw_names.close();

        System.out.println("X");
    }

    private static void format_uneven_ensgs(File f_output_rnaseq_sorted, File f_output_rnaseq_formatted) throws IOException {

        HashSet<String> available_ENSGs = new HashSet<>();

        HashMap<File,HashSet<String>> files_avail_ENSGs = new HashMap<>();
        HashMap<String,HashSet<File>> avail_ENSGs_files = new HashMap<>();

        int count_files = 0;


        for(File f_sorted : f_output_rnaseq_sorted.listFiles())
        {
            if(f_sorted.isDirectory())
            {
                for(File f_file : f_sorted.listFiles())
                {
                    if(f_file.isFile())
                    {
                        count_files++;

                        HashSet<String> avail_ensgs_this_file = new HashSet<>();
                        files_avail_ENSGs.put(f_file,avail_ensgs_this_file);

                        BufferedReader br = new BufferedReader(new FileReader(f_file));
                        String line= "";
                        while((line=br.readLine())!=null)
                        {
                            String[] split = line.split("\t");

                            String ensg = split[0].split("\\.")[0];
                            if(!ensg.startsWith("ENSG"))
                                continue;
                            available_ENSGs.add(ensg);
                            avail_ensgs_this_file.add(ensg);

                            if(avail_ENSGs_files.containsKey(ensg))
                            {
                                avail_ENSGs_files.get(ensg).add(f_file);
                            }
                            else
                            {
                                HashSet<File> avail_file = new HashSet<>();
                                avail_file.add(f_file);

                                avail_ENSGs_files.put(ensg,avail_file);
                            }

                        }
                        br.close();
                    }
                }
            }
        }

        ArrayList<String> removable_ensgs = new ArrayList<>();

        for(String k_ensg : avail_ENSGs_files.keySet())
        {
            if(k_ensg.startsWith("__"))
            {
                removable_ensgs.add(k_ensg);
                continue;
            }

            if(avail_ENSGs_files.get(k_ensg).size()!=88)
            {
                System.out.println("SIZE: " + k_ensg + " - " + avail_ENSGs_files.get(k_ensg).size());

                //removable_ensgs.add(k_ensg);
            }

            if(avail_ENSGs_files.get(k_ensg).size()==count_files)
            {
                System.out.println("in all files: " + k_ensg);
            }

            /*
            if(avail_ENSGs_files.get(k_ensg).size()!=count_files)
            {
                removable_ensgs.add(k_ensg);
            }*/
        }

        for(String k_ensg : removable_ensgs)
        {
            avail_ENSGs_files.remove(k_ensg);
        }

        ArrayList<String> ensgs_ordered = new ArrayList<>(avail_ENSGs_files.keySet());
        Collections.sort(ensgs_ordered);

        for(String key_ensg: ensgs_ordered)
        {
            if(!key_ensg.startsWith("ENSG"))
            {
                System.out.println("FAIL!!");
            }
        }

        HashMap<String,Integer> ensg_to_place = new HashMap<>();
        for(int i = 0; i < ensgs_ordered.size(); i++)
        {
            ensg_to_place.put(ensgs_ordered.get(i),i);
        }

        for(File f_sorted : f_output_rnaseq_sorted.listFiles())
        {
            if (f_sorted.isDirectory())
            {
                File f_output = new File(f_output_rnaseq_formatted.getAbsolutePath()+File.separator+f_sorted.getName());
                f_output.mkdirs();

                for (File f_file : f_sorted.listFiles())
                {
                    if (f_file.isFile())
                    {
                        HashMap<String,Boolean> already_got_ENSG = new HashMap<>();
                        HashMap<Integer,String> position_to_count = new HashMap<>();

                        for(String key_ensg : avail_ENSGs_files.keySet())
                        {
                            already_got_ENSG.put(key_ensg,false);
                        }

                        String[] rep_name_split= f_file.getName().split("_");
                        String replicate=rep_name_split[rep_name_split.length-1].split("\\.")[0];

                        String[] split_sorted = f_sorted.getName().split("_");
                        String groupname = "";
                        for(String k_split_sorted:split_sorted)
                        {
                            groupname+=k_split_sorted.substring(0, 1).toUpperCase() + k_split_sorted.substring(1);
                        }

                        String rep_name= groupname+"_Rep"+replicate+".txt";

                        File f_output_file = new File(f_output.getAbsolutePath()+File.separator+rep_name);

                        BufferedReader br = new BufferedReader(new FileReader(f_file));

                        String line = "";
                        while((line=br.readLine())!=null)
                        {
                            String[] split = line.split("\t");
                            String ensg = split[0].split("\\.")[0];
                            if(already_got_ENSG.containsKey(ensg))
                            {
                                if(already_got_ENSG.get(ensg))
                                {
                                    System.out.println("DOUBLE " + ensg);
                                    continue;
                                }
                                else
                                {
                                    already_got_ENSG.put(ensg,true);

                                    int place_in_file = ensg_to_place.get(ensg);

                                    position_to_count.put(place_in_file,split[1]);
                                }
                            }
                        }
                        br.close();

                        BufferedWriter bw = new BufferedWriter(new FileWriter(f_output_file));
                        bw.write(rep_name);
                        bw.newLine();

                        for(int i = 0 ; i < ensgs_ordered.size(); i++)
                        {
                            if(!position_to_count.containsKey(i))
                            {
                                String ensg = ensgs_ordered.get(i);

                                System.out.println("X");
                            }
                            if(!position_to_count.containsKey(i))
                            {
                                System.out.println("X");
                            }
                            bw.write(position_to_count.get(i));
                            bw.newLine();
                        }

                        bw.close();
                    }
                }
            }
        }

        File f_output_names_formatted = new File(f_output_rnaseq_formatted.getAbsolutePath()+File.separator+"names_formatted.txt");
        BufferedWriter bw_names = new BufferedWriter(new FileWriter(f_output_names_formatted));
        bw_names.write("Geneid");
        bw_names.newLine();
        for(int i = 0; i < ensgs_ordered.size();i++)
        {
            bw_names.write(ensgs_ordered.get(i));
            bw_names.newLine();
        }
        bw_names.close();

        System.out.println("X");
    }

    private static void sort(File f_input_root, File f_output_sorted_root, File f_output_rnaseq, HashSet<String> groups, String other_samples, File f_input_clinical, File f_input_gdc_sample_sheet, File f_output_mapping, HashMap<String,HashSet<String>> ident_folder,File f_output_rnaseq_sorted, File f_output_rnaseq_formatted, File f_output_pca) throws Exception
    {

        if(!ident_folder.isEmpty() && !execute_option.equals(""))
        {
            //create new folder for old data
            //move old data into new folder
            //create old file structure

            File f_output_old = new File(f_output_sorted_root.getAbsolutePath()+File.separator+"RNA_SEQ"+File.separator+"old_before_"+execute_option);
            f_output_old.mkdir();

            String command_move_sorted ="mv " + f_output_rnaseq_sorted + " " + f_output_old.getAbsolutePath();
            String command_move_formatted="mv " + f_output_rnaseq_formatted + " " + f_output_old.getAbsolutePath();
            String command_move_pca="mv " + f_output_pca + " " + f_output_old.getAbsolutePath();
            String command_move_mapping="mv " + f_output_mapping + " " + f_output_old.getAbsolutePath();

            Process child_move_sorted = Runtime.getRuntime().exec(command_move_sorted);
            int code_move_sorted = child_move_sorted.waitFor();
            switch (code_move_sorted){
                case 0:
                    break;
                case 1:
                    String message = child_move_sorted.getErrorStream().toString();
                    throw new Exception(message);
            }

            Process child_move_formatted = Runtime.getRuntime().exec(command_move_formatted);
            int code_move_formatted = child_move_formatted.waitFor();
            switch (code_move_formatted){
                case 0:
                    break;
                case 1:
                    String message = child_move_formatted.getErrorStream().toString();
                    throw new Exception(message);
            }

            Process child_move_pca = Runtime.getRuntime().exec(command_move_pca);
            int code_move_pca = child_move_pca.waitFor();
            switch (code_move_pca){
                case 0:
                    break;
                case 1:
                    String message = child_move_pca.getErrorStream().toString();
                    throw new Exception(message);
            }

            Process child_move_mapping = Runtime.getRuntime().exec(command_move_mapping);
            int code_move_mapping = child_move_mapping.waitFor();
            switch (code_move_mapping){
                case 0:
                    break;
                case 1:
                    String message = child_move_mapping.getErrorStream().toString();
                    throw new Exception(message);
            }

            f_output_rnaseq_sorted.mkdir();
            f_output_rnaseq_formatted.mkdir();
            f_output_pca.mkdir();
        }

        //create output files
        for(String k_group : groups)
        {
            File f_ouput = new File(f_output_rnaseq.getAbsolutePath()+File.separator+k_group);
            f_ouput.mkdirs();
        }
        File f_output_others = new File(f_output_rnaseq.getAbsolutePath()+File.separator+other_samples);
        f_output_others.mkdirs();

        HashMap<String,String> file_id_name = new HashMap<>();
        HashMap<String,String> file_name_id = new HashMap<>();
        HashMap<String,String> file_sample_id = new HashMap<>();
        HashMap<String,String> file_sample_name = new HashMap<>();
        HashMap<String,String> file_id_sample = new HashMap<>();
        HashMap<String,String> file_name_sample = new HashMap<>();
        HashMap<String,String> sample_subtype = new HashMap<>();
        HashMap<String,String> name_subtype = new HashMap<>();
        HashMap<String,String> sample_name_caseID = new HashMap<>();
        HashMap<String,String> file_id_sample_id = new HashMap<>();

        BufferedReader br_gdc_sample_sheet = new BufferedReader(new FileReader(f_input_gdc_sample_sheet));
        String line_gdc_sample_sheet = br_gdc_sample_sheet.readLine();
        String[] header_gdc_sample_sheet = line_gdc_sample_sheet.split("\t");
        while((line_gdc_sample_sheet= br_gdc_sample_sheet.readLine())!=null)
        {
            String[] split = line_gdc_sample_sheet.split("\t");
            String file_id = split[0];
            String file_name=split[1];
            String sample_name = split[5];
            String sample_id = split[6];

            file_id_sample_id.put(file_id,sample_id);
            file_id_name.put(file_id,file_name);
            file_name_id.put(file_name,file_id);
            file_sample_id.put(sample_name,file_id);
            file_sample_name.put(sample_name,file_name);
            file_id_sample.put(file_id,sample_name);
            file_name_sample.put(file_name,sample_name);

        }
        br_gdc_sample_sheet.close();

        BufferedReader br_clinical = new BufferedReader(new FileReader(f_input_clinical));
        String line_clinical = br_clinical.readLine();
        String[] header_clinical = line_clinical.split("\t");
        while((line_clinical= br_clinical.readLine())!=null)
        {
            String[] split = line_clinical.split("\t");
            String subtype = split[107];
            String caseID = split[0];
            String sample_name = split[1];

            sample_subtype.put(sample_name,subtype);
            name_subtype.put(caseID,subtype);
            sample_name_caseID.put(sample_name,caseID);
        }
        br_clinical.close();

        StringBuilder sb_mappings = new StringBuilder();
        sb_mappings.append("sample_name\tcase_submitter_id\tcase_ID\tgdc_fileID\tgdc_file_name\tsubtype\tspecific_subtype\tsample_id\n");

        for(File f_file_id: f_input_root.listFiles())
        {
            if (f_file_id.isDirectory())
            {
                boolean execute_this_file = true;

                String id = f_file_id.getName();


                if(!execute_option.equals(""))
                {
                    if(!ident_folder.isEmpty())
                    {
                        for(String key_ident : ident_folder.keySet())
                        {
                            if(execute_option.equals(key_ident))
                            {
                                continue;
                            }

                            if(ident_folder.get(key_ident).contains(id))
                            {
                                execute_this_file=false;
                                break;
                            }
                        }
                    }
                }

                if(!execute_this_file)
                {
                    continue;
                }

                String name = file_id_name.get(id);
                String sample = file_id_sample.get(id);
                String sample_id = file_id_sample_id.get(id);


                if(!sample_subtype.containsKey(sample))
                {
                    continue;
                }
                String subtype = sample_subtype.get(sample);
                String caseID = sample_name_caseID.get(sample);

                String subtype_name = "";
                String subtype_formatted ="";

                for(String k_group : groups)
                {
                    if(subtype.matches(".*"+k_group+".*"))
                    {
                        subtype_name=subtype;
                        subtype_formatted=k_group;
                    }
                }

                boolean other_than_group = false;

                if(subtype_name.equals(""))
                {
                    subtype_name=subtype;
                    subtype_formatted=other_samples;
                    other_than_group=true;
                }

                for(File f_intern: f_file_id.listFiles())
                {
                    if(f_intern.isFile())
                    {

                        String command = "cp " + f_intern.getAbsolutePath() + " ";
                        String command_move = "mv ";
                        String move_name="";

                        if(other_than_group)
                        {
                            File f_output_dir = new File(f_output_rnaseq.getAbsolutePath()+File.separator+subtype_formatted);

                            int sample_number = f_output_dir.listFiles().length;

                            move_name = other_samples + "_sample_"+sample_number + ".counts";

                            command+= f_output_dir.getAbsolutePath();
                            command_move += f_output_dir.getAbsolutePath()+File.separator+f_intern.getName()+" " + f_output_dir.getAbsolutePath()+File.separator+move_name;
                        }
                        else
                        {
                            File f_output_dir = new File(f_output_rnaseq.getAbsolutePath()+File.separator+subtype_formatted);

                            int sample_number = f_output_dir.listFiles().length;

                            move_name = subtype_formatted + "_sample_"+sample_number + ".counts";

                            command+= f_output_dir.getAbsolutePath();
                            command_move += f_output_dir.getAbsolutePath()+File.separator+f_intern.getName()+" " + f_output_dir.getAbsolutePath()+File.separator+move_name;
                        }

                        sb_mappings.append(move_name);
                        sb_mappings.append("\t");
                        sb_mappings.append(sample);
                        sb_mappings.append("\t");
                        sb_mappings.append(id);
                        sb_mappings.append("\t");
                        sb_mappings.append(caseID);
                        sb_mappings.append("\t");
                        sb_mappings.append(file_id_name.get(id));
                        sb_mappings.append("\t");
                        sb_mappings.append(subtype_formatted);
                        sb_mappings.append("\t");
                        sb_mappings.append(subtype_name);
                        sb_mappings.append("\t");
                        sb_mappings.append(sample_id);
                        sb_mappings.append("\n");

                        //copy file
                        Process child = Runtime.getRuntime().exec(command);
                        int code = child.waitFor();
                        switch (code){
                            case 0:
                                break;
                            case 1:
                                String message = child.getErrorStream().toString();
                                throw new Exception(message);
                        }

                        //rename file
                        Process child_move = Runtime.getRuntime().exec(command_move);
                        int code_move = child_move.waitFor();
                        switch (code_move){
                            case 0:
                                break;
                            case 1:
                                String message = child_move.getErrorStream().toString();
                                throw new Exception(message);
                        }
                    }
                }
            }
        }


        BufferedWriter bw_mapping = new BufferedWriter(new FileWriter(f_output_mapping));
        bw_mapping.write(sb_mappings.toString());
        bw_mapping.close();


        System.out.println("X");
    }

    private static boolean unpack(File input_root) throws Exception {

        boolean correct = false;

        String command_edited = "find "+input_root.getAbsolutePath()+" -name \\*.gz -exec sh -c 'echo {}; gzip -d {}' \\;";

        File f_script = new File(input_root.getParentFile().getAbsolutePath()+File.separator+"unpack.sh");
        BufferedWriter bw = new BufferedWriter(new FileWriter(f_script));
        bw.write("#!/bin/bash");
        bw.newLine();
        bw.write(command_edited);
        bw.close();

        //System.out.println("EXECUTE THIS COMMAND BY YOUR OWN");
        //System.out.println(command_edited);

        String command_shell ="sh " + f_script.getAbsolutePath();

        Process child = Runtime.getRuntime().exec(command_shell);
        int code = child.waitFor();
        switch (code){
            case 0:
                correct=true;
                break;
            case 1:
                String message = child.getErrorStream().toString();
                correct=false;
                //throw new Exception(message);
        }

        return correct;
    }
}
