package com2pose;

import util.Options_intern;

import javax.swing.*;
import java.io.*;
import java.util.*;

public class COM2POSE_lib
{
    Options_intern options_intern;

    public COM2POSE_lib(Options_intern options_intern)
    {
        this.options_intern = options_intern;
    }

    public void run_and_postprocess_DESeq2()
    {

    }

    /**
     *create DESeq2 scripts based on input directory for DESeq2 - each group against each group, save intermediate steps and R Scripts
     */
    public void create_DESeq2_scripts() throws IOException {

        HashMap<Integer,String> row_ensg_name=new HashMap<>();

        BufferedReader br_ensg_per_line= new BufferedReader(new FileReader(new File(options_intern.deseq2_input_gene_id)));
        String line_ensg_per_line="";
        line_ensg_per_line=br_ensg_per_line.readLine();
        int count_ensg_lines = 0;
        while((line_ensg_per_line=br_ensg_per_line.readLine())!=null)
        {
            row_ensg_name.put(count_ensg_lines,line_ensg_per_line);
            count_ensg_lines++;
        }
        br_ensg_per_line.close();

        File output_intermediate_steps = new File(options_intern.com2pose_working_directory+File.separator+"DESeq2_preprocessing");
        if(output_intermediate_steps.exists())
        {
            System.err.println("Working directory was already used - please us another one or empty this one completely");
            //TODO: after debugging use system exit !!!
            //System.exit(1);
        }
        output_intermediate_steps.mkdir();
        File output_inter_steps_combined=new File(options_intern.com2pose_working_directory+File.separator+"DESeq2_preprocessing"+File.separator+"combined");
        output_inter_steps_combined.mkdir();
        File output_inter_steps_single=new File(options_intern.com2pose_working_directory+File.separator+"DESeq2_preprocessing"+File.separator+"single");
        output_inter_steps_single.mkdir();
        File folder = new File(options_intern.deseq2_input_directory);

        HashMap<String, HashSet<String>> timepoints_samples= new HashMap<>();

        //CREATE SINGLES
        for(File fileDir:folder.listFiles())
        {
            if(fileDir.isDirectory())
            {
                File group = new File(output_inter_steps_single.getAbsolutePath()+File.separator+fileDir.getName());
                group.mkdir();

                HashSet<String> x = new HashSet<>();
                timepoints_samples.put(group.getName(),x);
                for(File sample : fileDir.listFiles())
                {
                    String name= sample.getName();
                    timepoints_samples.get(group.getName()).add(name);

                    BufferedReader br = new BufferedReader(new FileReader(sample));
                    BufferedWriter bw = new BufferedWriter(new FileWriter(group.getAbsolutePath()+File.separator+name));
                    String line = br.readLine();
                    bw.write(line);
                    bw.newLine();
                    int count=0;
                    while((line=br.readLine())!=null)
                    {
                        bw.write(row_ensg_name.get(count)+"\t"+line);
                        bw.newLine();
                        count++;
                    }
                    br.close();
                    bw.close();
                }
            }
        }
        //CREATE ALL COMBINED FILES
        HashSet<String> already_combined=new HashSet<>();
        for(String k: timepoints_samples.keySet())
        {
            for(String kk: timepoints_samples.keySet())
            {
                String key1=k+"_"+kk;
                String key2=kk+"_"+k;

                if(already_combined.contains(key1)||already_combined.contains(key2)||k.equals(kk))
                {
                    continue;
                }

                already_combined.add(key1);

                HashSet<String> samples_group1 = timepoints_samples.get(k);
                HashSet<String> samples_group2 = timepoints_samples.get(kk);

                File combined_file = new File(output_inter_steps_combined.getAbsolutePath()+File.separator+key1);
                combined_file.mkdir();

                StringBuilder sb_header = new StringBuilder();
                sb_header.append("geneID");

                ArrayList<BufferedReader> bufferedReaders = new ArrayList<>();
                for(String s: samples_group1)
                {
                    sb_header.append("\t");
                    sb_header.append(s);
                    BufferedReader br = new BufferedReader(new FileReader(new File(output_inter_steps_single+File.separator+k+File.separator+s)));
                    br.readLine();
                    bufferedReaders.add(br);
                }

                for(String s: samples_group2)
                {
                    sb_header.append("\t");
                    sb_header.append(s);
                    BufferedReader br = new BufferedReader(new FileReader(new File(output_inter_steps_single+File.separator+kk+File.separator+s)));
                    br.readLine();
                    bufferedReaders.add(br);
                }

                BufferedWriter bw = new BufferedWriter(new FileWriter(new File(output_inter_steps_combined+File.separator+key1+File.separator+key1+".csv")));
                bw.write(sb_header.toString());
                bw.newLine();

                BufferedReader br_ensg = new BufferedReader(new FileReader(new File(options_intern.deseq2_input_gene_id)));
                String line = "";
                br_ensg.readLine();
                while((line=br_ensg.readLine())!=null)
                {
                    StringBuilder sb = new StringBuilder();
                    sb.append(line);
                    for(BufferedReader br : bufferedReaders)
                    {
                        String line_intern = br.readLine();

                        String[] split = line_intern.split("\t");
                        sb.append("\t");
                        sb.append(split[1]);
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br_ensg.close();
                bw.close();
            }
        }

        //CREATE RScripts for every folder in here
        File r_scripts = new File(options_intern.com2pose_working_directory+File.separator+"DESeq2_R_scripts");
        r_scripts.mkdir();
        for(File dir : output_inter_steps_combined.listFiles())
        {
            if(dir.isDirectory())
            {
                String name_rscript = dir.getName()+".R";

                String group1 = dir.getName().split("_")[0];
                String group2 = dir.getName().split("_")[1];


                BufferedWriter bw = new BufferedWriter(new FileWriter(new File(r_scripts.getAbsolutePath()+File.separator+name_rscript)));

                StringBuilder sb = new StringBuilder();
                sb.append("require(\"DESeq2\")\n");
                sb.append("library(DESeq2)\n");
                sb.append("#"+group1+" VS "+group2+"\n");

                BufferedReader br_header = new BufferedReader(new FileReader(dir.getAbsolutePath()+File.separator+dir.getName()+".csv"));
                String line_header = br_header.readLine();
                br_header.close();

                ArrayList<String> groups = new ArrayList<>();
                ArrayList<String> samples = new ArrayList<>();

                String[] split_header= line_header.split("\t");
                for(int i = 1; i < split_header.length;i++)
                {
                    String[] split_samples = split_header[i].split("-");
                    String[] split_groups = split_samples[0].split("_");

                    samples.add(split_samples[0]);
                    groups.add(split_groups[0]);
                }

                sb.append("metadata_df<-data.frame(sample_id = c(");
                int count_s=0;
                for(String s: samples)
                {
                    if(count_s>0)
                    {
                        sb.append(" ,\"");
                    }
                    else
                    {
                        sb.append("\"");
                    }
                    sb.append(s);
                    sb.append("\"");
                    count_s++;
                }
                sb.append("), group = c(");
                count_s=0;
                for(String s: groups)
                {
                    if(count_s>0)
                    {
                        sb.append(" ,\"");
                    }
                    else
                    {
                        sb.append("\"");
                    }
                    sb.append(s);
                    sb.append("\"");
                    count_s++;
                }
                sb.append("))\n");

                sb.append("input_groups = \"");
                sb.append(dir.getName());
                sb.append("\"\n");

                sb.append("group_one = strsplit(input_groups, \"_\")[[1]][1]\n" +
                        "group_two = strsplit(input_groups, \"_\")[[1]][2]\n");

                sb.append("rownames(metadata_df) <- metadata_df$sample_id\n" +
                        "metadata_df$sample_id <- NULL\n");

                sb.append("count_path = \" ");
                sb.append(dir.getAbsolutePath());
                sb.append(File.separator+dir.getName());
                sb.append(".csv\"\n");
                sb.append("count_df = read.csv(count_path, sep = \"\\t\", header = T, row.names = 1)\n");

                sb.append("dds <- DESeqDataSetFromMatrix(countData=count_df, \n");
                sb.append("                              colData=metadata_df, \n");
                sb.append("                              design=~group)\n");

                if(options_intern.deseq2_count_threshold>0)
                {
                    sb.append("threshold = ");
                    sb.append(options_intern.deseq2_count_threshold);
                    sb.append("\n");
                    sb.append("keep <- rowSums(counts(dds)) >= threshold\n");
                    sb.append("dds <- dds[keep,]\n");
                }

                File output_deseq2 = new File(options_intern.com2pose_working_directory+File.separator+"DESeq2_output_raw");
                output_deseq2.mkdir();

                sb.append("output_path = \"");
                sb.append(output_deseq2.getAbsolutePath());
                sb.append(File.separator);
                sb.append(dir.getName());
                sb.append("_raw_output_deseq2.tsv");
                sb.append("\"\n");

                sb.append("dds <- DESeq(dds)\n");
                sb.append("res <- results(dds)\n");
                sb.append("summary(res)\n");
                sb.append("write.table(res,file=output_path,sep=\"\\t\")\n");


                bw.write(sb.toString());
                bw.close();
            }
        }

        System.out.println("X");


    }

    /**
     * read config file
     * @throws IOException
     */
    public void read_config_file() throws IOException {

        BufferedReader br = new BufferedReader(new FileReader(new File(options_intern.config_data_path)));
        String line= "";
        while((line=br.readLine())!=null)
        {
            if(line.startsWith("#"))
            {
                continue;
            }

            String[] split=line.split("=");

            switch (split[0])
            {
                case "com2pose_working_directory":
                    options_intern.com2pose_working_directory=split[1].substring(1,split[1].length()-1);
                    break;
                case "deseq2_input_directory":
                    options_intern.deseq2_input_directory=split[1].substring(1,split[1].length()-1);
                    break;
                case "deseq2_input_gene_id":
                    options_intern.deseq2_input_gene_id=split[1].substring(1,split[1].length()-1);
                    break;
                case "deseq2_count_threshold":
                    options_intern.deseq2_count_threshold=Integer.parseInt(split[1]);
                    break;
                case "tepic_input_directory":
                    options_intern.tepic_input_directory=split[1].substring(1,split[1].length()-1);
                    break;
                case "tepic_input_ref_genome":
                    options_intern.tepic_input_ref_genome=split[1].substring(1,split[1].length()-1);
                    break;
                case "tepic_path_pwms":
                    options_intern.tepic_path_pwms=split[1].substring(1,split[1].length()-1);
                    break;
                case "tepic_cores":
                    options_intern.tepic_cores=Integer.parseInt(split[1]);
                    break;
                case "tepic_bed_chr_sign":
                    options_intern.tepic_bed_chr_sign=split[1].substring(1,split[1].length()-1);
                    break;
                case "tepic_column_bedfile":
                    options_intern.tepic_column_bedfile=Integer.parseInt(split[1]);
                    break;
                case "tepic_gene_annot":
                    options_intern.tepic_gene_annot=split[1].substring(1,split[1].length()-1);
                    break;
                case "tepic_window_size":
                    options_intern.tepic_window_size=Integer.parseInt(split[1]);
                    break;
                case "tepic_onlyDNasePeaks":
                    options_intern.tepic_onlyDNasePeaks=split[1].substring(1,split[1].length()-1);
                    break;
                case "tepic_exponential_decay":
                    options_intern.tepic_exponential_decay=Boolean.parseBoolean(split[1]);
                    break;
                case "tepic_not_norm_peak_length":
                    options_intern.tepic_not_norm_peak_length=Boolean.parseBoolean(split[1]);
                    break;
                case "tepic_not_generated":
                    options_intern.tepic_not_generated=Boolean.parseBoolean(split[1]);
                    break;
                case "tepic_original_decay":
                    options_intern.tepic_original_decay=Boolean.parseBoolean(split[1]);
                    break;
                case "tepic_psems_length":
                    options_intern.tepic_psems_length=split[1].substring(1,split[1].length()-1);
                    break;
                case "tepic_entire_gene_body":
                    options_intern.tepic_entire_gene_body=Boolean.parseBoolean(split[1]);
                    break;
                case "tepic_zipped":
                    options_intern.tepic_zipped=Boolean.parseBoolean(split[1]);
                    break;
                case "tepic_background_seq":
                    options_intern.tepic_background_seq=split[1].substring(1,split[1].length()-1);
                    break;
                case "tepic_2bit":
                    options_intern.tepic_2bit=split[1].substring(1,split[1].length()-1);
                    break;
                case "tepic_pvalue":
                    options_intern.tepic_pvalue=Double.parseDouble(split[1]);
                    break;
                case "tepic_minutes_per_chr":
                    options_intern.tepic_minutes_per_chr=Integer.parseInt(split[1]);
                    break;
                case "tepic_chr_prefix":
                    options_intern.tepic_chr_prefix=Boolean.parseBoolean(split[1]);
                    break;
                case "tepic_transcript_based":
                    options_intern.tepic_transcript_based=Boolean.parseBoolean(split[1]);
                    break;
                case "tepic_loop_list":
                    options_intern.tepic_loop_list=split[1].substring(1,split[1].length()-1);
                    break;
                case "tepic_loop_windows":
                    options_intern.tepic_loop_windows=Integer.parseInt(split[1]);
                    break;
                case "tepic_only_peak_features":
                    options_intern.tepic_only_peak_features=Boolean.parseBoolean(split[1]);
                    break;
                case "tepic_tpm_cutoff":
                    options_intern.tepic_tpm_cutoff=Integer.parseInt(split[1]);
                    break;
                case "tepic_ensg_symbol":
                    options_intern.tepic_ensg_symbol=split[1].substring(1,split[1].length()-1);
                    break;
                default:
                    System.err.println("Misformed cfg file - please use template of: /COM2POSE/config_templates/com2pose_template.cfg");
                    System.err.println("Do not delete unused parameters in config data!");
                    System.exit(1);
            }
        }
        br.close();


    }


    /*
    public void run()
    {


        //process peak file or if directory all files in directory
        File folder = new File(options_intern.input_dir_peaks);
        run_all_files_subdirs(folder);

    }


    private void run_single(File file_in)
    {
        //The main program

    }


    private void run_all_files_subdirs(File folder)
    {
        if(folder.isDirectory())
        {
            for(File fileDir : folder.listFiles())
            {
                if(fileDir.isDirectory()&& options_intern.run_all_subdirectories)
                {
                    run_all_files_subdirs(fileDir);
                }
                else if (!fileDir.isDirectory())
                {
                    run_single(fileDir);
                }
            }
        }
        else
        {
            run_single(folder);
        }

    }*/



}
