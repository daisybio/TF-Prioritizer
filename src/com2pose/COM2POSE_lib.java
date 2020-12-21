package com2pose;
import util.Logger;
import util.Options_intern;
import java.io.*;
import java.util.*;

public class COM2POSE_lib
{
    Options_intern options_intern;
    Logger logger;

    //HashMap<String,HashMap<String,HashSet<String>>> tepic_groups_histoneModifications_samples = new HashMap<>();

    public COM2POSE_lib(Options_intern options_intern) throws IOException {
        this.options_intern = options_intern;

        logger = new Logger(options_intern.write_to_logfile, options_intern.com2pose_working_directory);

        logger.logLine("#########################################");
        logger.logLine("############## COM2POSE #################");
        logger.logLine("#########################################");
        logger.logLine("Working directory set to: " + options_intern.com2pose_working_directory);
        logger.logLine("COM2POSE path set to: "+ options_intern.path_to_COM2POSE);
    }

    public void run_tepic() throws Exception {
        logger.logLine("Start TEPIC.sh");

        String command = "bash";
        String tepic_path = " " + options_intern.path_to_COM2POSE+File.separator+"ext"+File.separator+"TEPIC"+File.separator+"TEPIC"+File.separator+"Code"+File.separator+"TEPIC.sh";
        command += tepic_path;
        command += " -g "+ options_intern.tepic_input_ref_genome;
        command += " -p "+ options_intern.tepic_path_pwms;


        String command_tail = "";
        if(options_intern.tepic_cores>1)
        {
            command_tail += " -c "+ options_intern.tepic_cores;
        }
        if(!options_intern.tepic_bed_chr_sign.equals(""))
        {
            command_tail += " -d "+ options_intern.tepic_bed_chr_sign;
        }
        if(options_intern.tepic_column_bedfile!=-1)
        {
            command_tail += " -n "+options_intern.tepic_column_bedfile;
        }
        if(!options_intern.tepic_gene_annot.equals(""))
        {
            command_tail += " -a " + options_intern.tepic_gene_annot;
        }
        if(options_intern.tepic_window_size!=50000)
        {
            command_tail += " -w " + options_intern.tepic_window_size;

        }
        if(!options_intern.tepic_onlyDNasePeaks.equals(""))
        {
            command_tail += " -f " + options_intern.tepic_onlyDNasePeaks;
        }
        if(options_intern.tepic_exponential_decay)
        {
            command_tail += " -e TRUE";
        }
        if(options_intern.tepic_not_norm_peak_length)
        {
            command_tail += " -l TRUE";
        }
        if(options_intern.tepic_not_generated)
        {
            command_tail += " -u TRUE";
        }
        if(options_intern.tepic_original_decay)
        {
            command_tail += " -x TRUE";
        }
        if(!options_intern.tepic_psems_length.equals(""))
        {
            command_tail += " -m "+ options_intern.tepic_psems_length;
        }
        if(options_intern.tepic_entire_gene_body)
        {
            command_tail += " -y TRUE";
        }
        if(options_intern.tepic_zipped)
        {
            command_tail += " -z TRUE";
        }
        if(!options_intern.tepic_2bit.equals(""))
        {
            command_tail += " -r " + options_intern.tepic_2bit;
        }
        if(options_intern.tepic_pvalue!=0.05)
        {
            command_tail += " -v " + options_intern.tepic_pvalue;

        }
        if(options_intern.tepic_minutes_per_chr!=3)
        {
            command_tail += " -i " + options_intern.tepic_minutes_per_chr;
        }
        if(options_intern.tepic_chr_prefix)
        {
            command_tail += " -j TRUE";
        }
        if(options_intern.tepic_transcript_based)
        {
            command_tail += " -t TRUE";
        }
        if(!options_intern.tepic_loop_list.equals(""))
        {
            command_tail += " -h " + options_intern.tepic_loop_list;
        }
        if(options_intern.tepic_loop_windows!=5000)
        {
            command_tail += " -s " + options_intern.tepic_loop_windows;
        }
        if(options_intern.tepic_only_peak_features)
        {
            command_tail+= " -q TRUE";
        }
        if(options_intern.tepic_tpm_cutoff>0)
        {
            command_tail += " -T " + options_intern.tepic_tpm_cutoff;
            command_tail += " -E " + options_intern.tepic_ensg_symbol;
            command_tail += " -A " + options_intern.deseq2_input_gene_id;
        }

        File output_TEPIC = new File(options_intern.com2pose_working_directory+File.separator+"TEPIC_output_raw");
        output_TEPIC.mkdir();


        File folder = new File(options_intern.tepic_input_directory);
        for(File dirGroup :folder.listFiles())
        {
            if(dirGroup.isDirectory())
            {
                File output_TEPIC_group = new File(output_TEPIC.getAbsolutePath()+File.separator+dirGroup.getName());
                output_TEPIC_group.mkdir();

                logger.logLine("[TEPIC] Start group "+ dirGroup.getName());
                for(File dirHM : dirGroup.listFiles())
                {
                    File output_TEPIC_group_hm = new File(output_TEPIC_group.getAbsolutePath()+File.separator+dirHM.getName());
                    output_TEPIC_group_hm.mkdir();

                    logger.logLine("[TEPIC] Start histone modification: "+ dirHM.getName());
                    for(File sample : dirHM.listFiles())
                    {
                        logger.logLine("[TEPIC] Start sample " + sample.getName());

                        File output_sample = new File(output_TEPIC_group_hm.getAbsolutePath()+File.separator+sample.getName());
                        output_sample.mkdir();

                        String command_sample = new String(command);

                        String output_dir_combined = output_sample.getAbsolutePath()+File.separator+output_sample.getName();

                        command_sample += " -b " + sample.getAbsolutePath();
                        command_sample += " -o " + output_dir_combined;

                        String command_tail_sample = new String(command_tail);
                        if(options_intern.tepic_tpm_cutoff>0)
                        {
                            String n_dir = options_intern.com2pose_working_directory+File.separator+"DESeq2_preprocessing"+File.separator+"single"+File.separator+dirGroup.getName()+"_meanCounts.txt";
                            command_tail_sample += " -G " + n_dir;
                        }

                        String command_execute = command_sample + command_tail_sample;
                        logger.logLine("[TEPIC] execute TEPIC with command line: " + command_execute);
                        Process child = Runtime.getRuntime().exec(command_execute);
                        int code = child.waitFor();
                        switch (code){
                            case 0:
                                break;
                            case 1:
                                String message = child.getErrorStream().toString();
                                throw new Exception(message);
                        }


                        System.out.println("X");

                    }
                }
            }
        }

        logger.logLine("Finished TEPIC.sh");
    }

    /**
     *run created DESeq2 scripts and postprocess for input into DYNAMITE
     */
    public void run_and_postprocess_DESeq2() throws Exception {

        logger.logLine("Start running DESeq2 RScripts");

        File folder = new File(options_intern.com2pose_working_directory+File.separator+"DESeq2_R_scripts");

        for (File dir:folder.listFiles())
        {
            if(!dir.isDirectory())
            {
                String command = "Rscript " + dir.getAbsolutePath();
                Process child = Runtime.getRuntime().exec(command);
                logger.logLine("[DESEQ2] Running script " + dir.getName());
                int code = child.waitFor();
                switch (code){
                    case 0:
                        break;
                    case 1:
                        String message = child.getErrorStream().toString();
                        throw new Exception(message);
                }
            }
        }

        logger.logLine("Finished running DESeq2 RScripts");
        logger.logLine("Start postprocessing DESeq2 data for input to DYNAMITE");


        File folder_results = new File(options_intern.com2pose_working_directory+File.separator+"DESeq2_output_raw");
        File output_file = new File(options_intern.com2pose_working_directory+File.separator+"DESeq2_output");
        output_file.mkdir();

        for(File res:folder_results.listFiles())
        {
            if(!res.isDirectory())
            {
                String name_res = res.getName();
                String[] split_name_res = name_res.split("_");
                String name_out = split_name_res[0]+"_"+split_name_res[1]+"_DYNAMITE.tsv";

                BufferedReader br = new BufferedReader(new FileReader(res));
                BufferedWriter bw = new BufferedWriter(new FileWriter(output_file.getAbsolutePath()+File.separator+name_out));

                bw.write("geneID\tlog2fc");
                bw.newLine();

                String line = br.readLine();

                while((line= br.readLine())!=null)
                {
                    String[] split = line.split("\t");
                    if(!split[2].equals("NA"))
                    {
                        StringBuilder sb = new StringBuilder();
                        sb.append(split[0].substring(1,split[0].length()-1));
                        sb.append("\t");
                        sb.append(split[2]);

                        bw.write(sb.toString());
                        bw.newLine();
                    }
                }
                bw.close();
                br.close();




            }
        }



        logger.logLine("Finished Start postprocessing DESeq2 data for input to DYNAMITE");

    }

    /**
     *create DESeq2 scripts based on input directory for DESeq2 - each group against each group, save intermediate steps and R Scripts
     */
    public void create_DESeq2_scripts() throws IOException {

        logger.logLine("Start preprocessing nfcore RNA-seq for DESeq2 input");
        HashMap<Integer,String> row_ensg_name=new HashMap<>();

        File ensg_names = new File(options_intern.deseq2_input_gene_id);
        BufferedReader br_ensg_per_line= new BufferedReader(new FileReader(ensg_names));
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
            logger.logLine("Working directory was already used - please us another one or empty this one completely");
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
                HashMap<Integer,Integer> mean_line_counts = new HashMap<>();
                int count_samples = 0;
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
                        int count_line = Integer.parseInt(line);
                        if(mean_line_counts.containsKey(count))
                        {
                            int z = mean_line_counts.get(count);
                            z+= count_line;
                            mean_line_counts.put(count,z);
                        }
                        else
                        {
                            mean_line_counts.put(count,count_line);
                        }
                        count++;
                    }
                    if(count!=count_ensg_lines)
                    {
                        logger.logLine("Error in nfcore RNA-seq data: File "+ sample.getName()+ " has not the same number of rows as in File "+ ensg_names.getName());
                        System.exit(1);
                    }
                    br.close();
                    bw.close();

                    count_samples++;
                }

                BufferedWriter bw_means = new BufferedWriter(new FileWriter(new File(output_inter_steps_single.getAbsolutePath()+File.separator+group.getName()+"_meanCounts.txt")));
                bw_means.write(group.getName()+"_MEANS");
                bw_means.newLine();
                for(int i = 0; i < count_ensg_lines; i++)
                {
                    int mean_count = mean_line_counts.get(i)/count_samples;
                    bw_means.write(""+mean_count);
                    bw_means.newLine();

                }
                bw_means.close();


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

        logger.logLine("Finished preprocessing of nfcore RNA-seq data - no errors detected");
        logger.logLine("Started creating RScripts for running DESeq2");

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
                sb.append("if (!requireNamespace(\"BiocManager\", quietly = TRUE))\n");
                sb.append("  install.packages(\"BiocManager\")\n");
                sb.append("if (!requireNamespace(\"DESeq2\", quietly = TRUE))\n");
                sb.append("  install.packages(\"DESeq2\")\n");
                sb.append("library(DESeq2)\n");
                sb.append("#"+group1+" VS "+group2+"\n");
                logger.logLine("[DESEQ2] Group " + group1 + " VS " + group2);

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

                sb.append("group_one = strsplit(input_groups, \"_\")[[1]][1]\n");
                sb.append("group_two = strsplit(input_groups, \"_\")[[1]][2]\n");

                sb.append("rownames(metadata_df) <- metadata_df$sample_id\n");
                sb.append("metadata_df$sample_id <- NULL\n");

                sb.append("count_path = \"");
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

        logger.logLine("Finished creating RScripts for running DESeq2");


    }

    /**
     * read config file
     */
    public void read_config_file() throws IOException {

        logger.logLine("Start reading config file at "+ options_intern.config_data_path);

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
                    options_intern.tepic_tpm_cutoff=Double.parseDouble(split[1]);
                    break;
                case "tepic_ensg_symbol":
                    options_intern.tepic_ensg_symbol=split[1].substring(1,split[1].length()-1);
                    break;
                default:
                    logger.logLine("Misformed cfg file - please use template of: /COM2POSE/config_templates/com2pose_template.cfg");
                    logger.logLine("Do not delete unused parameters in config data!");
                    System.exit(1);
            }
        }
        br.close();

        boolean all_set = checkOptions();
        logger.logLine("Check config file parameters for validity");
        if(!all_set)
        {
            logger.logLine("Not all REQuired options set. Please set them in config file");
            logger.logLine("Aborting COM2POSE");
            System.exit(1);
        }
        else
        {
            logger.logLine("Parameters in config file valid");
        }

        logger.logLine("Reading config file finished - no errors detected");


    }


    private boolean checkOptions() throws IOException {

        boolean all_set = true;

        /**
         * DESEQ2 options
         */

        if(options_intern.deseq2_input_directory.equals(""))
        {
            logger.logLine("[DESEQ2] input directory is not given");
            all_set=false;
        }
        else
        {
            File f = new File(options_intern.deseq2_input_directory);
            if(!f.exists())
            {
                logger.logLine("[DESEQ2] input path does not exist!");
                all_set=false;
            }
        }

        if(options_intern.deseq2_input_gene_id.equals(""))
        {
            logger.logLine("[DESEQ2] gene ID file from nfcore RNA-seq is not given");
            all_set=false;
        }
        else
        {
            File f = new File(options_intern.deseq2_input_gene_id);
            if(!f.exists())
            {
                logger.logLine("[DESEQ2] gene ID file from nfcore RNA-seq path does not exist!");
                all_set=false;
            }
        }


        /**
         * TEPIC options
         */
        if(options_intern.tepic_input_ref_genome.equals(""))
        {
            logger.logLine("[TEPIC] Reference genome path is not given");
            all_set=false;
        }
        else
        {
            File f = new File(options_intern.tepic_input_ref_genome);
            if(!f.exists())
            {
                logger.logLine("[TEPIC] Reference genome path does not exist!");
                all_set=false;
            }
        }
        if(options_intern.tepic_input_directory.equals(""))
        {
            logger.logLine("[TEPIC] nfcore ChIP-seq data directory is not given");
            all_set=false;
        }
        else
        {
            File f = new File(options_intern.tepic_input_directory);
            if(!f.exists())
            {
                logger.logLine("[TEPIC] nfcore ChIP-seq data path does not exist!");
                all_set=false;
            }
        }

        if(options_intern.tepic_path_pwms.equals(""))
        {
            logger.logLine("[TEPIC] position specific energy matrix directory is not given");
            all_set=false;
        }
        else
        {
            File f = new File(options_intern.tepic_path_pwms);
            if(!f.exists())
            {
                logger.logLine("[TEPIC] position specific energy matrix path does not exist!");
                all_set=false;
            }
        }

        if(options_intern.tepic_tpm_cutoff>0)
        {
            //check for gene annotation
            if(options_intern.tepic_gene_annot.equals(""))
            {
                logger.logLine("[TEPIC] TPM cutoff set, but no annotation file is given");
                all_set=false;
            }
            else
            {
                File f = new File(options_intern.tepic_gene_annot);
                if(!f.exists())
                {
                    logger.logLine("[TEPIC] TPM cutoff set and gene annotation file path does not exist!");
                    all_set=false;
                }
            }

            //check for map ensg symbol
            if(options_intern.tepic_ensg_symbol.equals(""))
            {
                logger.logLine("[TEPIC] TPM cutoff set, but no map of ENSG to Gene Symbol is given!");
            }
            else
            {
                File f = new File(options_intern.tepic_ensg_symbol);
                if(!f.exists())
                {
                    logger.logLine("[TEPIC] TPM cutoff set and map of ENSG to Gene Symbol file path does not exist!");
                    all_set=false;
                }
            }
        }

        if(!options_intern.tepic_background_seq.equals("") && !options_intern.tepic_2bit.equals(""))
        {
            logger.logLine("[TEPIC] parameters: tepic_background_seq and tepic_2bit are mutually exclusive");
            all_set=false;
        }

        if(options_intern.tepic_column_bedfile!=-1 && !options_intern.tepic_bed_chr_sign.equals(""))
        {
            logger.logLine("[TEPIC] parameters: tepic_column_bedfile and tepic_bed_chr_sign are mutually exclusive");
            all_set=false;
        }

        if(!options_intern.tepic_bed_chr_sign.equals(""))
        {
            File f = new File(options_intern.tepic_bed_chr_sign);
            if(!f.exists())
            {
                logger.logLine("[TEPIC] tepic_bed_chr_sign file path does not exists!");
                all_set=false;
            }
        }

        if(!options_intern.tepic_psems_length.equals(""))
        {
            File f = new File(options_intern.tepic_psems_length);
            if(!f.exists())
            {
                logger.logLine("[TEPIC] tepic_psems_length file path does not exists!");
                all_set=false;
            }
        }

        if(!options_intern.tepic_loop_list.equals(""))
        {
            File f = new File(options_intern.tepic_loop_list);
            if(!f.exists())
            {
                logger.logLine("[TEPIC] tepic_loop_list file path does not exists!");
                all_set=false;
            }
        }



        /**
         * DYNAMITE OPTIONS
         */


        return all_set;
    }




}
