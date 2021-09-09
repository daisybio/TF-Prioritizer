package tcga_atac_seq_helpers;

import java.io.*;
import java.util.HashMap;
import java.util.HashSet;

public class atac_seq_download_convert_sort_identify_remove {

    public static void main(String[] args) throws Exception {

        //COUNT WITH BASH DISTINCT VALUES OF A COLUMN
        //awk -F '\t' '{print $108}' clinical.tsv | sort | uniq -c | sort -nr

        String input_gdc_client = "/home/markus/Downloads/gdc-client_v1.6.1_Ubuntu_x64";

        /**
         * TGCT block
         */
        atac_seq_download_cancer_blocks tgct = new atac_seq_download_cancer_blocks();
        {
            tgct.name="TGCT";
            tgct.input_manifest="/home/markus/Downloads/gdc_manifest_20210716_083509.txt";
            tgct.input_gdc_sample_sheet="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/TGCT/METADATA/gdc_sample_sheet.2021-07-16.tsv";
            tgct.input_clinical="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/TGCT/METADATA/clinical.cart.2021-07-16/clinical.tsv";
            tgct.output_directory="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/TGCT";
            tgct.groups = new HashSet<>();
            tgct.groups.add("Seminoma");
            tgct.other_samples = "NonSeminoma";

            //execute_all(tgct);
        }


        /**
         * BLCA block
         */
        atac_seq_download_cancer_blocks blca = new atac_seq_download_cancer_blocks();
        {
            blca.name="BLCA";
            blca.input_manifest="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/03_BLCA/META_DATA/ATAC_SEQ/gdc_manifest_20210720_141715.txt";
            blca.input_clinical="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/03_BLCA/META_DATA/ATAC_SEQ/clinical.cart.2021-07-20/clinical.tsv";
            blca.input_gdc_sample_sheet="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/03_BLCA/META_DATA/ATAC_SEQ/gdc_sample_sheet.2021-07-20.tsv";
            blca.output_directory="/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/03_BLCA/ATAC_SEQ";
            blca.groups = new HashSet<>();
            blca.groups.add("Papillary");
            blca.other_samples = "NonPapillary";
            //execute_all(blca);
        }

        /**
         * 03 BRCA BLOCK
         */
        atac_seq_download_cancer_blocks brca = new atac_seq_download_cancer_blocks();
        {
            //COUNT WITH BASH DISTINCT VALUES OF A COLUMN
            //awk -F '\t' '{print $108}' clinical.tsv | sort | uniq -c | sort -nr

            brca.name="03_BRCA";
            brca.input_manifest="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/03_BRCA/METADATA/ATAC_SEQ/gdc_manifest_20210721_125249.txt";
            brca.input_clinical="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/03_BRCA/METADATA/ATAC_SEQ/clinical.tsv";
            brca.input_gdc_sample_sheet="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/03_BRCA/METADATA/ATAC_SEQ/gdc_sample_sheet.2021-07-30.tsv";
            brca.output_directory="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/03_BRCA/ATAC_SEQ";
            brca.groups = new HashSet<>();
            brca.groups.add("Infiltrating");
            brca.groups.add("Lobular");
            brca.groups.add("Metaplastic");
            brca.other_samples = "Others";
            execute_all(brca);
        }

        /**
         * 04 COAD BLOCK
         */
        atac_seq_download_cancer_blocks coad = new atac_seq_download_cancer_blocks();
        {
            //COUNT WITH BASH DISTINCT VALUES OF A COLUMN
            //awk -F '\t' '{print $108}' clinical.tsv | sort | uniq -c | sort -nr

            coad.name="04_COAD";
            coad.input_manifest="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/04_COAD/METADATA/ATAC_SEQ/gdc_manifest_20210722_083841.txt";
            coad.input_clinical="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/04_COAD/METADATA/ATAC_SEQ/clinical.cart.2021-07-22/clinical.tsv";
            coad.input_gdc_sample_sheet="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/04_COAD/METADATA/ATAC_SEQ/gdc_sample_sheet.2021-07-22.tsv";
            coad.output_directory="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/04_COAD/ATAC_SEQ";
            coad.groups = new HashSet<>();
            coad.groups.add("Adenocarcinoma");
            coad.other_samples = "NonAdenocarcinoma";
            execute_all(coad);
        }

        /**
         * 05 KIRP BLOCK
         */
        atac_seq_download_cancer_blocks kirp = new atac_seq_download_cancer_blocks();
        {
            //COUNT WITH BASH DISTINCT VALUES OF A COLUMN
            //awk -F '\t' '{print $108}' clinical.tsv | sort | uniq -c | sort -nr

            kirp.name="05_KIRP";
            kirp.input_manifest="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/05_KIRP/METADATA/ATAC_SEQ/gdc_manifest_20210722_084501.txt";
            kirp.input_clinical="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/05_KIRP/METADATA/ATAC_SEQ/clinical.cart.2021-07-22/clinical.tsv";
            kirp.input_gdc_sample_sheet="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/05_KIRP/METADATA/ATAC_SEQ/gdc_sample_sheet.2021-07-22.tsv";
            kirp.output_directory="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/05_KIRP/ATAC_SEQ";
            kirp.groups = new HashSet<>();
            kirp.groups.add("Papillary");
            kirp.other_samples = "NonPapillary";
            execute_all(kirp);
        }

        /**
         * 06_PRAD BLOCK
         */
        atac_seq_download_cancer_blocks prad = new atac_seq_download_cancer_blocks();
        {
            //COUNT WITH BASH DISTINCT VALUES OF A COLUMN
            //awk -F '\t' '{print $108}' clinical.tsv | sort | uniq -c | sort -nr

            prad.name="06_PRAD";
            prad.input_manifest="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/06_PRAD/METADATA/ATAC_SEQ/gdc_manifest_20210722_084614.txt";
            prad.input_clinical="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/06_PRAD/METADATA/ATAC_SEQ/clinical.cart.2021-07-22/clinical.tsv";
            prad.input_gdc_sample_sheet="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/06_PRAD/METADATA/ATAC_SEQ/gdc_sample_sheet.2021-07-22.tsv";
            prad.output_directory="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/06_PRAD/ATAC_SEQ";
            prad.groups = new HashSet<>();
            prad.groups.add("Adenocarcinoma");
            prad.other_samples = "Infiltrating";
            execute_all(prad);
        }

        /**
         * 07_LUAD BLOCK
         */
        atac_seq_download_cancer_blocks luad = new atac_seq_download_cancer_blocks();
        {
            //COUNT WITH BASH DISTINCT VALUES OF A COLUMN
            //awk -F '\t' '{print $108}' clinical.tsv | sort | uniq -c | sort -nr

            luad.name="07_LUAD";
            luad.input_manifest="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/07_LUAD/METADATA/ATAC_SEQ/gdc_manifest_20210722_084723.txt";
            luad.input_clinical="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/07_LUAD/METADATA/ATAC_SEQ/clinical.cart.2021-07-22/clinical.tsv";
            luad.input_gdc_sample_sheet="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/07_LUAD/METADATA/ATAC_SEQ/gdc_sample_sheet.2021-07-22.tsv";
            luad.output_directory="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/07_LUAD/ATAC_SEQ";
            luad.groups = new HashSet<>();
            luad.groups.add("Adenocarcinoma");
            luad.groups.add("Acinar");
            luad.groups.add("Mucinous");
            luad.other_samples = "PapillaryAdenocarcinoma";
            execute_all(luad);
        }

        /**
         * 08_STAD BLOCK
         */
        atac_seq_download_cancer_blocks stad = new atac_seq_download_cancer_blocks();
        {
            //COUNT WITH BASH DISTINCT VALUES OF A COLUMN
            //awk -F '\t' '{print $108}' clinical.tsv | sort | uniq -c | sort -nr

            stad.name="08_STAD";
            stad.input_manifest="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/08_STAD/METADATA/ATAC_SEQ/gdc_manifest_20210722_084854.txt";
            stad.input_clinical="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/08_STAD/METADATA/ATAC_SEQ/clinical.cart.2021-07-22/clinical.tsv";
            stad.input_gdc_sample_sheet="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/08_STAD/METADATA/ATAC_SEQ/gdc_sample_sheet.2021-07-22.tsv";
            stad.output_directory="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/08_STAD/ATAC_SEQ";
            stad.groups = new HashSet<>();
            stad.groups.add("Adenocarcinoma");
            stad.groups.add("Mucinous");
            stad.groups.add("Tubular");
            stad.other_samples = "NonAdenocarcinoma";
            execute_all(stad);
        }

        /**
         * 09_ESCA
         */
        atac_seq_download_cancer_blocks esca = new atac_seq_download_cancer_blocks();
        {
            //COUNT WITH BASH DISTINCT VALUES OF A COLUMN
            //awk -F '\t' '{print $108}' clinical.tsv | sort | uniq -c | sort -nr

            esca.name="09_ESCA";
            esca.input_manifest="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/09_ESCA/METADATA/ATAC_SEQ/gdc_manifest_20210722_084954.txt";
            esca.input_clinical="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/09_ESCA/METADATA/ATAC_SEQ/clinical.cart.2021-07-22/clinical.tsv";
            esca.input_gdc_sample_sheet="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/09_ESCA/METADATA/ATAC_SEQ/gdc_sample_sheet.2021-07-22.tsv";
            esca.output_directory="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/09_ESCA/ATAC_SEQ";
            esca.groups = new HashSet<>();
            esca.groups.add("Squamous");
            esca.other_samples = "Adenocarcinoma";
            execute_all(esca);
        }

        /**
         * 10_LIHC
         */
        atac_seq_download_cancer_blocks lihc = new atac_seq_download_cancer_blocks();
        {
            //COUNT WITH BASH DISTINCT VALUES OF A COLUMN
            //awk -F '\t' '{print $108}' clinical.tsv | sort | uniq -c | sort -nr

            lihc.name="10_LIHC";
            lihc.input_manifest="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/10_LIHC/METADATA/ATAC_SEQ/gdc_manifest_20210722_085059.txt";
            lihc.input_clinical="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/10_LIHC/METADATA/ATAC_SEQ/clinical.cart.2021-07-22/clinical.tsv";
            lihc.input_gdc_sample_sheet="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/10_LIHC/METADATA/ATAC_SEQ/gdc_sample_sheet.2021-07-22.tsv";
            lihc.output_directory="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/10_LIHC/ATAC_SEQ";
            lihc.groups = new HashSet<>();
            lihc.groups.add("Hepatocellular");
            lihc.other_samples = "NonHepatocellular";
            execute_all(lihc);
        }

        /**
         * 11_KIRC
         */
        atac_seq_download_cancer_blocks kirc = new atac_seq_download_cancer_blocks();
        {
            //COUNT WITH BASH DISTINCT VALUES OF A COLUMN
            //awk -F '\t' '{print $108}' clinical.tsv | sort | uniq -c | sort -nr

            kirc.name="11_KIRC";
            kirc.input_manifest="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/11_KIRC/METADATA/ATAC_SEQ/gdc_manifest_20210722_085201.txt";
            kirc.input_clinical="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/11_KIRC/METADATA/ATAC_SEQ/clinical.cart.2021-07-22/clinical.tsv";
            kirc.input_gdc_sample_sheet="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/11_KIRC/METADATA/ATAC_SEQ/gdc_sample_sheet.2021-07-22.tsv";
            kirc.output_directory="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/11_KIRC/ATAC_SEQ";
            kirc.groups = new HashSet<>();
            kirc.groups.add("Clear");
            kirc.other_samples = "Renal,";
            execute_all(kirc);
        }

        /**
         * 12_LUSC
         */
        atac_seq_download_cancer_blocks lusc = new atac_seq_download_cancer_blocks();
        {
            //COUNT WITH BASH DISTINCT VALUES OF A COLUMN
            //awk -F '\t' '{print $108}' clinical.tsv | sort | uniq -c | sort -nr

            lusc.name="12_LUSC";
            lusc.input_manifest="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/12_LUSC/METADATA/ATAC_SEQ/gdc_manifest_20210722_085303.txt";
            lusc.input_clinical="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/12_LUSC/METADATA/ATAC_SEQ/clinical.cart.2021-07-22/clinical.tsv";
            lusc.input_gdc_sample_sheet="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/12_LUSC/METADATA/ATAC_SEQ/gdc_sample_sheet.2021-07-22.tsv";
            lusc.output_directory="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/12_LUSC/ATAC_SEQ";
            lusc.groups = new HashSet<>();
            lusc.groups.add("keratinizing");
            lusc.other_samples = "Squamous";
            execute_all(lusc);
        }

        /**
         * 13_THCA
         */
        atac_seq_download_cancer_blocks thca = new atac_seq_download_cancer_blocks();
        {
            //COUNT WITH BASH DISTINCT VALUES OF A COLUMN
            //awk -F '\t' '{print $108}' clinical.tsv | sort | uniq -c | sort -nr

            thca.name="13_THCA";
            thca.input_manifest="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/13_THCA/METADATA/ATAC_SEQ/gdc_manifest_20210722_085403.txt";
            thca.input_clinical="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/13_THCA/METADATA/ATAC_SEQ/clinical.cart.2021-07-22/clinical.tsv";
            thca.input_gdc_sample_sheet="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/13_THCA/METADATA/ATAC_SEQ/gdc_sample_sheet.2021-07-22.tsv";
            thca.output_directory="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/13_THCA/ATAC_SEQ";
            thca.groups = new HashSet<>();
            thca.groups.add("Papillary");
            thca.groups.add("Oxyphilic");
            thca.other_samples = "Follicular";
            execute_all(thca);
        }

        /**
         * 14_LGG
         */
        atac_seq_download_cancer_blocks lgg = new atac_seq_download_cancer_blocks();
        {
            //COUNT WITH BASH DISTINCT VALUES OF A COLUMN
            //awk -F '\t' '{print $108}' clinical.tsv | sort | uniq -c | sort -nr

            lgg.name="14_LGG";
            lgg.input_manifest="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/14_LGG/METADATA/ATAC_SEQ/gdc_manifest_20210722_085513.txt";
            lgg.input_clinical="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/14_LGG/METADATA/ATAC_SEQ/clinical.cart.2021-07-22/clinical.tsv";
            lgg.input_gdc_sample_sheet="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/14_LGG/METADATA/ATAC_SEQ/gdc_sample_sheet.2021-07-22.tsv";
            lgg.output_directory="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/14_LGG/ATAC_SEQ";
            lgg.groups = new HashSet<>();
            lgg.groups.add("Oligodendroglioma");
            lgg.groups.add("Astrocytoma");
            lgg.other_samples = "Mixed";
            execute_all(lgg);
        }

        /**
         * 15_SKCM
         */
        atac_seq_download_cancer_blocks skcm = new atac_seq_download_cancer_blocks();
        {
            //COUNT WITH BASH DISTINCT VALUES OF A COLUMN
            //awk -F '\t' '{print $108}' clinical.tsv | sort | uniq -c | sort -nr

            skcm.name="15_SKCM";
            skcm.input_manifest="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/15_SKCM/METADATA/ATAC_SEQ/gdc_manifest_20210722_085610.txt";
            skcm.input_clinical="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/15_SKCM/METADATA/ATAC_SEQ/clinical.cart.2021-07-22/clinical.tsv";
            skcm.input_gdc_sample_sheet="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/15_SKCM/METADATA/ATAC_SEQ/gdc_sample_sheet.2021-07-22.tsv";
            skcm.output_directory="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/15_SKCM/ATAC_SEQ";
            skcm.groups = new HashSet<>();
            skcm.groups.add("Malignant");
            skcm.other_samples = "Epithelioid";
            execute_all(skcm);
        }

        /**
         * 16_UCEC
         */
        atac_seq_download_cancer_blocks ucec = new atac_seq_download_cancer_blocks();
        {
            //COUNT WITH BASH DISTINCT VALUES OF A COLUMN
            //awk -F '\t' '{print $108}' clinical.tsv | sort | uniq -c | sort -nr

            ucec.name="16_UCEC";
            ucec.input_manifest="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/16_UCEC/METADATA/ATAC_SEQ/gdc_manifest_20210722_085753.txt";
            ucec.input_clinical="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/16_UCEC/METADATA/ATAC_SEQ/clinical.cart.2021-07-22/clinical.tsv";
            ucec.input_gdc_sample_sheet="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/16_UCEC/METADATA/ATAC_SEQ/gdc_sample_sheet.2021-07-22.tsv";
            ucec.output_directory="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/16_UCEC/ATAC_SEQ";
            ucec.groups = new HashSet<>();
            ucec.groups.add("Endometrioid");
            ucec.other_samples = "SerousCystadenocarcinoma";
            execute_all(ucec);
        }

        /**
         * 17_ACC
         */
        atac_seq_download_cancer_blocks acc = new atac_seq_download_cancer_blocks();
        {
            //COUNT WITH BASH DISTINCT VALUES OF A COLUMN
            //awk -F '\t' '{print $108}' clinical.tsv | sort | uniq -c | sort -nr

            acc.name="17_ACC";
            acc.input_manifest="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/17_ACC/METADATA/ATAC_SEQ/gdc_manifest_20210722_085849.txt";
            acc.input_clinical="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/17_ACC/METADATA/ATAC_SEQ/clinical.cart.2021-07-22/clinical.tsv";
            acc.input_gdc_sample_sheet="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/17_ACC/METADATA/ATAC_SEQ/gdc_sample_sheet.2021-07-22.tsv";
            acc.output_directory="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/17_ACC/ATAC_SEQ";
            acc.groups = new HashSet<>();
            acc.groups.add("Adrenal");
            acc.other_samples = "NonAdrenal";
            execute_all(acc);
        }

        /**
         * 18_GBM
         */
        atac_seq_download_cancer_blocks gbm = new atac_seq_download_cancer_blocks();
        {
            //COUNT WITH BASH DISTINCT VALUES OF A COLUMN
            //awk -F '\t' '{print $108}' clinical.tsv | sort | uniq -c | sort -nr

            gbm.name="18_GBM";
            gbm.input_manifest="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/18_GBM/METADATA/ATAC_SEQ/gdc_manifest_20210722_085953.txt";
            gbm.input_clinical="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/18_GBM/METADATA/ATAC_SEQ/clinical.cart.2021-07-22/clinical.tsv";
            gbm.input_gdc_sample_sheet="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/18_GBM/METADATA/ATAC_SEQ/gdc_sample_sheet.2021-07-22.tsv";
            gbm.output_directory="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/18_GBM/ATAC_SEQ";
            gbm.groups = new HashSet<>();
            gbm.groups.add("Glioblastoma");
            gbm.other_samples = "NonGlioblastoma";
            execute_all(gbm);
        }

        /**
         * 19_HNSC
         */
        atac_seq_download_cancer_blocks hnsc = new atac_seq_download_cancer_blocks();
        {
            //COUNT WITH BASH DISTINCT VALUES OF A COLUMN
            //awk -F '\t' '{print $108}' clinical.tsv | sort | uniq -c | sort -nr

            //NOTE: THIS ONE BY HAND non-kreatin. and kreatin.!!
            hnsc.name="19_HNSC";
            hnsc.input_manifest="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/19_HNSC/METADATA/ATAC_SEQ/gdc_manifest_20210722_090106.txt";
            hnsc.input_clinical="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/19_HNSC/METADATA/ATAC_SEQ/clinical.cart.2021-07-22/clinical.tsv";
            hnsc.input_gdc_sample_sheet="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/19_HNSC/METADATA/ATAC_SEQ/gdc_sample_sheet.2021-07-22.tsv";
            hnsc.output_directory="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/19_HNSC/ATAC_SEQ";
            hnsc.groups = new HashSet<>();
            hnsc.groups.add("Basaloid");
            hnsc.groups.add("keratinizing");
            hnsc.other_samples = "Squamous";
            execute_all(hnsc);
        }

        /**
         * 20_PCPG
         */
        atac_seq_download_cancer_blocks pcpg = new atac_seq_download_cancer_blocks();
        {
            //COUNT WITH BASH DISTINCT VALUES OF A COLUMN
            //awk -F '\t' '{print $108}' clinical.tsv | sort | uniq -c | sort -nr

            pcpg.name="20_PCPG";
            pcpg.input_manifest="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/20_PCPG/METADATA/ATAC_SEQ/gdc_manifest_20210722_090212.txt";
            pcpg.input_clinical="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/20_PCPG/METADATA/ATAC_SEQ/clinical.cart.2021-07-22/clinical.tsv";
            pcpg.input_gdc_sample_sheet="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/20_PCPG/METADATA/ATAC_SEQ/gdc_sample_sheet.2021-07-22.tsv";
            pcpg.output_directory="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/20_PCPG/ATAC_SEQ";
            pcpg.groups = new HashSet<>();
            pcpg.groups.add("malignant");
            pcpg.other_samples = "NonMalignant";
            execute_all(pcpg);
        }

        /**
         * 21_MESO
         */
        atac_seq_download_cancer_blocks meso = new atac_seq_download_cancer_blocks();
        {
            //COUNT WITH BASH DISTINCT VALUES OF A COLUMN
            //awk -F '\t' '{print $108}' clinical.tsv | sort | uniq -c | sort -nr

            meso.name="21_MESO";
            meso.input_manifest="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/21_MESO/METADATA/ATAC_SEQ/gdc_manifest_20210722_090309.txt";
            meso.input_clinical="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/21_MESO/METADATA/ATAC_SEQ/clinical.cart.2021-07-22/clinical.tsv";
            meso.input_gdc_sample_sheet="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/21_MESO/METADATA/ATAC_SEQ/gdc_sample_sheet.2021-07-22.tsv";
            meso.output_directory="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/21_MESO/ATAC_SEQ";
            meso.groups = new HashSet<>();
            meso.groups.add("biphasic");
            meso.groups.add("Epithelioid");
            meso.other_samples = "Mesothelioma";
            execute_all(meso);
        }

        /**
         * 22_CHOL
         */
        atac_seq_download_cancer_blocks chol = new atac_seq_download_cancer_blocks();
        {
            //COUNT WITH BASH DISTINCT VALUES OF A COLUMN
            //awk -F '\t' '{print $108}' clinical.tsv | sort | uniq -c | sort -nr

            chol.name="22_CHOL";
            chol.input_manifest="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/22_CHOL/METADATA/ATAC_SEQ/gdc_manifest_20210722_090359.txt";
            chol.input_clinical="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/22_CHOL/METADATA/ATAC_SEQ/clinical.cart.2021-07-22/clinical.tsv";
            chol.input_gdc_sample_sheet="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/22_CHOL/METADATA/ATAC_SEQ/gdc_sample_sheet.2021-07-22.tsv";
            chol.output_directory="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/22_CHOL/ATAC_SEQ";
            chol.groups = new HashSet<>();
            chol.groups.add("Cholangiocarcinoma");
            chol.other_samples = "NonCholangiocarcinoma";
            execute_all(chol);
        }

        /**
         * 23_CESC
         */
        atac_seq_download_cancer_blocks cesc = new atac_seq_download_cancer_blocks();
        {
            //COUNT WITH BASH DISTINCT VALUES OF A COLUMN
            //awk -F '\t' '{print $108}' clinical.tsv | sort | uniq -c | sort -nr

            cesc.name="23_CESC";
            cesc.input_manifest="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/23_CESC/METADATA/ATAC_SEQ/gdc_manifest_20210722_090450.txt";
            cesc.input_clinical="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/23_CESC/METADATA/ATAC_SEQ/clinical.cart.2021-07-22/clinical.tsv";
            cesc.input_gdc_sample_sheet="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/23_CESC/METADATA/ATAC_SEQ/gdc_sample_sheet.2021-07-22.tsv";
            cesc.output_directory="/nfs/data/COM2POSE/TCGA/08_TCGA_ATAC_SEQ/23_CESC/ATAC_SEQ";
            cesc.groups = new HashSet<>();
            cesc.groups.add("Squamous");
            cesc.other_samples = "NonSquamous";
            execute_all(cesc);
        }

        /**
         * BLOCK
         */
        atac_seq_download_cancer_blocks example = new atac_seq_download_cancer_blocks();
        {
            //COUNT WITH BASH DISTINCT VALUES OF A COLUMN
            //awk -F '\t' '{print $108}' clinical.tsv | sort | uniq -c | sort -nr

            example.name="EXAMPLE";
            example.input_manifest="";
            example.input_clinical="";
            example.input_gdc_sample_sheet="";
            example.output_directory="/ATAC_SEQ";
            example.groups = new HashSet<>();
            example.groups.add("");
            example.other_samples = "";
            //execute_all(example);
        }

    }

    public static void execute_all(atac_seq_download_cancer_blocks block) throws Exception {

        System.out.println("EXECUTING DOWNLOADS FOR CANCER: " + block.name);

        String input_gdc_client= block.input_gdc_client;
        String input_manifest= block.input_manifest;
        String input_gdc_sample_sheet= block.input_gdc_sample_sheet;
        String input_clinical= block.input_clinical;
        String output_directory= block.output_directory;

        HashSet<String> groups = block.groups;
        String other_samples= block.other_samples;

        File f_output_directory = new File(output_directory);
        File f_output_directory_atac_seq = new File(f_output_directory.getAbsolutePath()+File.separator+"ATAC_SEQ");
        f_output_directory_atac_seq.mkdirs();
        File f_output_dir_raw_data = new File(f_output_directory_atac_seq.getAbsolutePath()+File.separator+"01_raw_data");
        f_output_dir_raw_data.mkdirs();
        File f_output_dir_sorted_data = new File(f_output_directory_atac_seq.getAbsolutePath()+File.separator+"02_sorted_data");
        f_output_dir_sorted_data.mkdirs();
        File f_output_dir_sorted_data_gapped = new File(f_output_directory_atac_seq.getAbsolutePath()+File.separator+"03_sorted_data_gapped");
        f_output_dir_sorted_data_gapped.mkdirs();
        File f_output_dir_sorted_data_rscript = new File(f_output_directory_atac_seq.getAbsolutePath()+File.separator+"04_sorted_data_rscript");
        f_output_dir_sorted_data_rscript.mkdirs();
        File f_output_dir_sorted_data_allpeaks= new File(f_output_directory_atac_seq.getAbsolutePath()+File.separator+"05_sorted_data_allpeaks");
        f_output_dir_sorted_data_allpeaks.mkdirs();

        //execute_download(input_manifest,input_gdc_client, f_output_dir_raw_data, block.input_token, block.name);

        //execute_peak_calling(f_output_dir_raw_data);

        File f_input_sample = new File(input_gdc_sample_sheet);
        if(!f_input_sample.exists())
        {
            System.out.println("CANNOT FIND SAMPLE SHEET");
            return;
        }
        File f_clinical = new File(input_clinical);
        if(!f_clinical.exists())
        {
            System.out.println("CANNOT FIND CLINICAL");
            return;
        }

        execute_sorting(input_gdc_sample_sheet,input_clinical,f_output_directory,f_output_dir_raw_data,f_output_dir_sorted_data,f_output_dir_sorted_data_gapped, groups, other_samples, f_output_dir_sorted_data_rscript, f_output_dir_sorted_data_allpeaks);

        //execute_chr_remove(f_output_dir_raw_data);

    }

    private static void execute_chr_remove(File f_output_dir_raw_data) {
        System.out.println("remove chr prefixes in .broadPeak files");
        String command = "find "+f_output_dir_raw_data.getAbsolutePath()+" -name \\*.broadPeak -exec sh -c 'echo {}; sed -i -e 's/chr//g' {}' \\;";
        System.out.println(command);
    }

    private static void execute_sorting(String input_gdc_sample_sheet, String input_clinical, File f_output_directory, File f_output_dir_raw_data, File f_output_dir_sorted_data, File f_output_dir_sorted_data_gapped, HashSet<String> groups, String other_samples, File f_output_dir_sorted_data_rscript, File f_output_dir_sorted_data_allpeaks) throws Exception {

        System.out.println("START SORTING");
        //create all output_dirs
        for(String k_groups : groups)
        {
            File f_out_sorted_data_parent = new File(f_output_dir_sorted_data.getAbsolutePath()+File.separator+k_groups);
            f_out_sorted_data_parent.mkdirs();

            File f_out_sorted_data= new File(f_out_sorted_data_parent.getAbsolutePath()+File.separator+"ATAC_SEQ");
            f_out_sorted_data.mkdirs();

            File f_out_gapped_parent = new File(f_output_dir_sorted_data_gapped.getAbsolutePath()+File.separator+k_groups);
            f_out_gapped_parent.mkdirs();

            File f_out_gapped = new File(f_out_gapped_parent.getAbsolutePath()+File.separator+"ATAC_SEQ");
            f_out_gapped.mkdirs();

            File f_out_rscript_parent = new File(f_output_dir_sorted_data_rscript.getAbsolutePath()+File.separator+k_groups);
            f_out_rscript_parent.mkdirs();

            File f_out_rscript= new File(f_out_rscript_parent.getAbsolutePath()+File.separator+"ATAC_SEQ");
            f_out_rscript.mkdirs();

            File f_out_allpeaks_parent = new File(f_output_dir_sorted_data_allpeaks.getAbsolutePath()+File.separator+k_groups);
            f_out_allpeaks_parent.mkdirs();

            File f_out_allpeaks = new File(f_out_allpeaks_parent.getAbsolutePath()+File.separator+"ATAC_SEQ");
            f_out_allpeaks.mkdirs();
        }

        File f_out_sorted_data_parent = new File(f_output_dir_sorted_data.getAbsolutePath()+File.separator+other_samples);
        f_out_sorted_data_parent.mkdirs();

        File f_out_sorted_data = new File(f_out_sorted_data_parent.getAbsolutePath()+File.separator+"ATAC_SEQ");
        f_out_sorted_data.mkdirs();

        File f_out_gapped_parent = new File(f_output_dir_sorted_data_gapped.getAbsolutePath()+File.separator+other_samples);
        f_out_gapped_parent.mkdirs();

        File f_out_gapped = new File(f_out_gapped_parent.getAbsolutePath()+File.separator+"ATAC_SEQ");
        f_out_gapped.mkdirs();

        File f_out_rscript_parent = new File(f_output_dir_sorted_data_rscript.getAbsolutePath()+File.separator+other_samples);
        f_out_rscript_parent.mkdirs();

        File f_out_rscript = new File(f_out_rscript_parent.getAbsolutePath()+File.separator+"ATAC_SEQ");
        f_out_rscript.mkdirs();

        File f_out_allpeaks_parent = new File(f_output_dir_sorted_data_allpeaks.getAbsolutePath()+File.separator+other_samples);
        f_out_allpeaks_parent.mkdirs();

        File f_out_allpeaks = new File(f_out_allpeaks_parent.getAbsolutePath()+File.separator+"ATAC_SEQ");
        f_out_allpeaks.mkdirs();


        File f_input_gdc_sample_sheet = new File(input_gdc_sample_sheet);
        File f_input_clinical = new File(input_clinical);

        HashMap<String,String> file_id_name = new HashMap<>();
        HashMap<String,String> file_name_id = new HashMap<>();
        HashMap<String,String> file_sample_id = new HashMap<>();
        HashMap<String,String> file_sample_name = new HashMap<>();
        HashMap<String,String> file_id_sample = new HashMap<>();
        HashMap<String,String> file_name_sample = new HashMap<>();
        HashMap<String,String> sample_subtype = new HashMap<>();
        HashMap<String,String> name_subtype = new HashMap<>();
        HashMap<String,String> sample_name_caseID = new HashMap<>();

        BufferedReader br_gdc_sample_sheet = new BufferedReader(new FileReader(f_input_gdc_sample_sheet));
        String line_gdc_sample_sheet = br_gdc_sample_sheet.readLine();
        String[] header_gdc_sample_sheet = line_gdc_sample_sheet.split("\t");
        while((line_gdc_sample_sheet= br_gdc_sample_sheet.readLine())!=null)
        {
            String[] split = line_gdc_sample_sheet.split("\t");
            String file_id = split[0];
            String file_name=split[1];
            String sample_name = split[5];

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
        sb_mappings.append("sample_name\tcase_submitter_id\tcase_ID\tgdc_fileID\tgdc_file_name\tsubtype\tspecific_subtype\n");

        for(File f_down_file : f_output_dir_raw_data.listFiles())
        {
            if(f_down_file.isDirectory())
            {
                String file_id = f_down_file.getName();
                if(file_id_sample.containsKey(file_id))
                {
                    String sample_name = file_id_sample.get(file_id);

                    String file_name_intern = file_id_name.get(file_id);
                    file_name_intern +="_file";

                    String subtype_name = sample_subtype.get(sample_name);
                    String subtype_sorted = "";

                    for(String key_subtypes : groups)
                    {
                        if(subtype_name.matches(".*"+key_subtypes+".*"))
                        {
                            subtype_sorted = key_subtypes;
                            subtype_name = key_subtypes;
                        }
                    }

                    boolean other_sample_here = false;
                    if(subtype_sorted.equals(""))
                    {
                        subtype_sorted=other_samples;
                        other_sample_here=true;
                    }

                    File f_macs2_files = new File(f_down_file.getAbsolutePath()+File.separator+file_name_intern);
                    if(f_macs2_files.exists())
                    {
                        if(f_macs2_files.isDirectory())
                        {
                            for(File f_macs2_outputs : f_macs2_files.listFiles())
                            {
                                String command = "cp " + f_macs2_outputs.getAbsolutePath() + " ";
                                String command_move = "mv ";

                                String name_f_move = f_macs2_outputs.getName();

                                if(name_f_move.matches(".*broadPeak.*"))
                                {
                                    String f_name_moved = "";

                                    if(other_sample_here)
                                    {
                                        int sample_number = f_out_sorted_data.listFiles().length;
                                        sample_number++;

                                        f_name_moved=other_samples+"_sample_"+sample_number+".broadPeak";
                                        command += f_out_sorted_data.getAbsolutePath();
                                        command_move+= f_out_sorted_data.getAbsolutePath()+File.separator+f_macs2_outputs.getName() + " " + f_out_sorted_data.getAbsolutePath()+File.separator+f_name_moved;

                                    }
                                    else
                                    {
                                        File f_output_group = new File(f_output_dir_sorted_data.getAbsolutePath()+File.separator+subtype_sorted+File.separator+"ATAC_SEQ");

                                        int sample_number = f_output_group.listFiles().length;
                                        sample_number++;

                                        f_name_moved=subtype_sorted+"_sample_"+sample_number+".broadPeak";
                                        command += f_output_group.getAbsolutePath();

                                        command_move+= f_output_group.getAbsolutePath()+File.separator+f_macs2_outputs.getName() + " " + f_output_group.getAbsolutePath()+File.separator+f_name_moved;

                                    }

                                    sb_mappings.append(f_name_moved);
                                    sb_mappings.append("\t");
                                    sb_mappings.append(sample_name);
                                    sb_mappings.append("\t");
                                    sb_mappings.append(file_id);
                                    sb_mappings.append("\t");
                                    sb_mappings.append(sample_name_caseID.get(sample_name));
                                    sb_mappings.append("\t");
                                    sb_mappings.append(file_id_name.get(file_id));
                                    sb_mappings.append("\t");
                                    sb_mappings.append(subtype_sorted);
                                    sb_mappings.append("\t");
                                    sb_mappings.append(subtype_name);
                                    sb_mappings.append("\n");


                                }

                                if(name_f_move.matches(".*gappedPeak.*"))
                                {
                                    String f_name_moved = "";

                                    if(other_sample_here)
                                    {
                                        int sample_number = f_out_gapped.listFiles().length;
                                        sample_number++;

                                        f_name_moved=other_samples+"_sample_"+sample_number+".gappedPeak";
                                        command += f_out_gapped.getAbsolutePath()+File.separator;
                                        command_move+= f_out_gapped.getAbsolutePath()+File.separator+f_macs2_outputs.getName() + " " + f_out_gapped.getAbsolutePath()+File.separator+f_name_moved;

                                    }
                                    else
                                    {
                                        File f_output_group = new File(f_output_dir_sorted_data_gapped.getAbsolutePath()+File.separator+subtype_sorted+File.separator+"ATAC_SEQ");

                                        int sample_number = f_output_group.listFiles().length;
                                        sample_number++;

                                        f_name_moved=subtype_sorted+"_sample_"+sample_number+".gappedPeak";
                                        command += f_output_group.getAbsolutePath();
                                        command_move+= f_output_group.getAbsolutePath()+File.separator+f_macs2_outputs.getName() + " " + f_output_group.getAbsolutePath()+File.separator+f_name_moved;

                                    }
                                }

                                if(name_f_move.matches(".*_model.*"))
                                {
                                    String f_name_moved = "";

                                    if(other_sample_here)
                                    {
                                        int sample_number = f_out_rscript.listFiles().length;
                                        sample_number++;

                                        f_name_moved=other_samples+"_sample_"+sample_number+"_model.R";
                                        command += f_out_rscript.getAbsolutePath();
                                        command_move+= f_out_rscript.getAbsolutePath()+File.separator+f_macs2_outputs.getName() + " " + f_out_rscript.getAbsolutePath()+File.separator+f_name_moved;


                                    }
                                    else
                                    {
                                        File f_output_group = new File(f_output_dir_sorted_data_rscript.getAbsolutePath()+File.separator+subtype_sorted+File.separator+"ATAC_SEQ");

                                        int sample_number = f_output_group.listFiles().length;
                                        sample_number++;

                                        f_name_moved=subtype_sorted+"_sample_"+sample_number+"_model.R";
                                        command += f_output_group.getAbsolutePath();
                                        command_move+= f_output_group.getAbsolutePath()+File.separator+f_macs2_outputs.getName() + " " + f_output_group.getAbsolutePath()+File.separator+f_name_moved;

                                    }
                                }

                                if(name_f_move.matches(".*_peaks.xls.*"))
                                {
                                    String f_name_moved = "";

                                    if(other_sample_here)
                                    {
                                        int sample_number = f_out_allpeaks.listFiles().length;
                                        sample_number++;

                                        f_name_moved=other_samples+"_sample_"+sample_number+"_peaks.xls";
                                        command += f_out_allpeaks.getAbsolutePath();
                                        command_move+= f_out_allpeaks.getAbsolutePath()+File.separator+f_macs2_outputs.getName() + " " + f_out_allpeaks.getAbsolutePath()+File.separator+f_name_moved;


                                    }
                                    else
                                    {
                                        File f_output_group = new File(f_output_dir_sorted_data_allpeaks.getAbsolutePath()+File.separator+subtype_sorted+File.separator+"ATAC_SEQ");

                                        int sample_number = f_output_group.listFiles().length;
                                        sample_number++;

                                        f_name_moved=subtype_sorted+"_sample_"+sample_number+"_peaks.xls";
                                        command += f_output_group.getAbsolutePath();
                                        command_move+= f_output_group.getAbsolutePath()+File.separator+f_macs2_outputs.getName() + " " + f_output_group.getAbsolutePath()+File.separator+f_name_moved;

                                    }
                                }

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

            }
        }

        File f_output_mapping = new File(f_output_directory.getAbsolutePath()+File.separator+"mapping_samples.csv");

        BufferedWriter bw_output_mapping = new BufferedWriter(new FileWriter(f_output_mapping));
        bw_output_mapping.append(sb_mappings.toString());
        bw_output_mapping.close();


        System.out.println("FINISHED SORTING");

    }

    private static void execute_peak_calling(File f_output_dir_raw_data) throws Exception
    {
        System.out.println("START PEAK CALLING");
        String command_mkdir = "find "+f_output_dir_raw_data.getAbsolutePath()+" -name \\*.bam -exec sh -c 'echo {}; mkdir {}_file' \\;";
        String command_edited = "find "+f_output_dir_raw_data.getAbsolutePath()+" -name \\*.bam -exec sh -c 'echo {}; macs2 callpeak --broad -t {} --outdir {}_file --format BAM' \\;";

        System.out.println(command_edited);
        /*
        Process child = Runtime.getRuntime().exec(command_edited);
        int code = child.waitFor();
        switch (code){
            case 0:
                break;
            case 1:
                String message = child.getErrorStream().toString();
                throw new Exception(message);
        }
        System.out.println("FINISHED PEAK CALLING");*/


    }

    private static void execute_download(String input_manifest, String input_gdc_client, File f_output_dir_raw_data, String input_token, String name) throws Exception
    {
        System.out.println("DOWNLOADING");
        File f_input_manifest = new File(input_manifest);

        String command_edited = input_gdc_client+ " download -m "+input_manifest+" --log-file "+f_input_manifest.getParentFile().getAbsolutePath()+File.separator+name+"_logfile.txt -t "+input_token+" -d " + f_output_dir_raw_data.getAbsolutePath();

        System.out.println(command_edited);
        /*
        Process child = Runtime.getRuntime().exec(command_edited);
        int code = child.waitFor();
        switch (code){
            case 0:
                break;
            case 1:
                String message = child.getErrorStream().toString();
                throw new Exception(message);
        }*/

        System.out.println("FINISHED DOWNLOADING");


    }
}
