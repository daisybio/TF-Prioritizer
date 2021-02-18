package com2pose;
import util.*;
import java.io.*;
import java.util.*;

public class COM2POSE_lib
{
    public Options_intern options_intern;
    public Logger logger;

    //only used in create_overview_website()!!
    ArrayList<File> threshold_folders = new ArrayList<>();
    boolean threshold_folders_filled = false;

    /**
     * constructor for analysis programms, you cannot run the pipeline here
     * @param options_intern
     * @param logger
     */
    public COM2POSE_lib(Options_intern options_intern, Logger logger)
    {
        this.options_intern = options_intern;
        this.logger=logger;
    }

    /**
     * constructor for the pipeline
     * @param options_intern
     * @throws IOException
     */
    public COM2POSE_lib(Options_intern options_intern) throws IOException {
        this.options_intern = options_intern;

        logger = new Logger(options_intern.write_to_logfile, options_intern.com2pose_working_directory);

        logger.logLine("#########################################");
        logger.logLine("############## COM2POSE #################");
        logger.logLine("#########################################");
        logger.logLine("Working directory set to: " + options_intern.com2pose_working_directory);
        logger.logLine("COM2POSE path set to: "+ options_intern.path_to_COM2POSE);
    }

    /**
     * CREATE HTML REPORT FOR DISTRIBUTION ANALYSIS
     */
    public void create_overview_html_report_distribution() throws IOException
    {
        logger.logLine("[DISTRIBUTION-ANALYSIS-HTML-REPORT] Create HTML report.");

        File f_output_website_root = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_out_website);

        File f_website_distr_analysis_html_folder = new File(f_output_website_root.getAbsolutePath()+File.separator+options_intern.folder_out_website_htmls_distribution_analysis);
        f_website_distr_analysis_html_folder.mkdir();
        File f_website_distr_analysis_html_folder_ALL = new File(f_website_distr_analysis_html_folder.getAbsolutePath()+File.separator+options_intern.folder_out_website_htmls_distribution_analysis_ALL);
        f_website_distr_analysis_html_folder_ALL.mkdir();
        File f_website_distr_analysis_html_folder_HM = new File(f_website_distr_analysis_html_folder.getAbsolutePath()+File.separator+options_intern.folder_out_website_htmls_distribution_analysis_HM);
        f_website_distr_analysis_html_folder_HM.mkdir();

        String html_header = get_header_html("HOME",options_intern.analysis_types_distribution_analysis);
        String html_tail = "</body>\n" +
                "</html>";

        File html_home_distribution_analysis = new File(f_output_website_root.getAbsolutePath()+File.separator+options_intern.html_report_home_regression_distribution_analysis);
        StringBuilder sb_home = new StringBuilder();

        sb_home.append(html_header);
        sb_home.append("<div class='w3-row-padding w3-padding-64 w3-container w3-content'><a href='' target='_blank'><button class='button_expandable'>ALL</button></a></div>");




        sb_home.append(html_tail);

        BufferedWriter bw_html_home = new BufferedWriter(new FileWriter(html_home_distribution_analysis));
        bw_html_home.write(sb_home.toString());
        bw_html_home.close();


        logger.logLine("[DISTRIBUTION-ANALYSIS-HTML-REPORT] Finished HTML report.");
    }

    /**
     * CREATE python scripts for plots, to only accept TFs which have a higher mean TF-TG-SCORE than the background distribution
     */
    public void create_distribution_plots() throws Exception {
        logger.logLine("[DISTRIBUTION-ANALYSIS-PLOTS] Start comparing background distribution to TF distribution and create plots for outstanding TFs");

        File f_analysis_distr_root = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_out_distribution);

        File f_analysis_distr_HM_level_input = new File(f_analysis_distr_root.getAbsolutePath()+File.separator+options_intern.folder_out_distribution_tf_tg_scores+File.separator+options_intern.folder_out_distribution_tf_tg_scores_background_distr+File.separator+options_intern.folder_out_distribution_tf_tg_scores_background_distr_HM);

        File f_root_plots_scripts = new File(f_analysis_distr_root.getAbsolutePath()+File.separator+options_intern.folder_out_distribution_plots_scripts);
        f_root_plots_scripts.mkdir();
        File f_root_plots_scripts_all = new File(f_root_plots_scripts.getAbsolutePath()+File.separator+options_intern.folder_out_distribution_plots_script_ALL);
        f_root_plots_scripts_all.mkdir();
        File f_root_plots_scripts_hm = new File(f_root_plots_scripts.getAbsolutePath()+File.separator+options_intern.folder_out_distribution_plots_scripts_HM);
        f_root_plots_scripts_hm.mkdir();


        File f_root_plots = new File(f_analysis_distr_root.getAbsolutePath()+File.separator+options_intern.folder_out_distribution_plots);
        f_root_plots.mkdir();
        File f_root_plots_all = new File(f_root_plots.getAbsolutePath()+File.separator+options_intern.folder_out_distribution_plots_ALL);
        f_root_plots_all.mkdir();
        File f_root_plots_hm = new File(f_root_plots.getAbsolutePath()+File.separator+options_intern.folder_out_distribution_plots_HM);
        f_root_plots_hm.mkdir();

        File f_scores_root = new File(f_analysis_distr_root.getAbsolutePath()+File.separator+options_intern.folder_out_distribution_tf_tg_scores);

        File f_background_distr_root = new File(f_scores_root.getAbsolutePath()+File.separator+options_intern.folder_out_distribution_tf_tg_scores_background_distr);
        File f_background_distr_tfs_ALL = new File(f_scores_root.getAbsolutePath()+File.separator+options_intern.folder_out_distribution_tf_tg_scores_tf_distributions+File.separator+options_intern.folder_out_distribution_tf_tg_scores_background_distr_ALL);
        File f_background_distr_tfs_HM = new File(f_scores_root.getAbsolutePath()+File.separator+options_intern.folder_out_distribution_tf_tg_scores_background_distr+File.separator+options_intern.folder_out_distribution_tf_tg_scores_background_distr_HM);

        File f_tf_distribution_ALL = new File(f_scores_root.getAbsolutePath()+File.separator+options_intern.folder_out_distribution_tf_tg_scores_tf_distributions+File.separator+options_intern.folder_out_distribution_tf_tg_scores_tf_distributions_ALL);
        File f_tf_distribution_HM = new File(f_scores_root.getAbsolutePath()+File.separator+options_intern.folder_out_distribution_tf_tg_scores_tf_distributions+File.separator+options_intern.folder_out_distribution_tf_tg_scores_tf_distributions_HM);

        File f_stats = new File(f_analysis_distr_root.getAbsolutePath()+File.separator+options_intern.folder_out_distribution_stats);
        f_stats.mkdir();
        File f_stats_ALL = new File(f_stats.getAbsolutePath()+File.separator+options_intern.folder_out_distribution_stats_ALL);
        f_stats_ALL.mkdir();
        File f_stats_HM = new File(f_stats.getAbsolutePath()+File.separator+options_intern.folder_out_distribution_stats_HM);
        f_stats_HM.mkdir();

        ArrayList<File> scripts_to_execute = new ArrayList<>();

        for(File fileDir_hm : f_analysis_distr_HM_level_input.listFiles())
        {
            if(fileDir_hm.isDirectory())
            {
                File f_output_script = new File(f_root_plots_scripts_hm.getAbsolutePath()+File.separator+fileDir_hm.getName());
                f_output_script.mkdir();
                File f_output_script_file = new File(f_output_script.getAbsolutePath()+File.separator+options_intern.file_suffix_distribution_analysis_python_script);

                File f_output_plot = new File(f_root_plots_hm.getAbsolutePath()+File.separator+fileDir_hm.getName());
                f_output_plot.mkdir();

                File f_output_stats = new File(f_stats_HM.getAbsolutePath()+File.separator+fileDir_hm.getName());
                f_output_stats.mkdir();

                File f_input_background_distr = new File(f_background_distr_tfs_HM.getAbsolutePath()+File.separator+fileDir_hm.getName()+File.separator+options_intern.file_suffix_distribution_analysis_distributions);
                File f_input_tf_distr = new File(f_tf_distribution_HM.getAbsolutePath()+File.separator+fileDir_hm.getName());

                scripts_to_execute.add(f_output_script_file);

                //File input_background_file, File input_tf_root, File output_plots, File output_script_file, File output_stats) throws IOException {

                write_python_script_distribution_analysis(f_input_background_distr,f_input_tf_distr,f_output_plot,f_output_script_file,f_output_stats);
            }
        }

        File all_background_distr= new File(f_background_distr_root.getAbsolutePath()+File.separator+options_intern.folder_out_distribution_tf_tg_scores_background_distr_ALL+File.separator+options_intern.file_suffix_distribution_analysis_distributions);
        File f_script_out_all = new File(f_root_plots_scripts_all.getAbsolutePath()+File.separator+options_intern.file_suffix_distribution_analysis_python_script);

        write_python_script_distribution_analysis(all_background_distr,f_tf_distribution_ALL,f_root_plots_all,f_script_out_all,f_stats_ALL);
        scripts_to_execute.add(f_script_out_all);

        for(File fileDir: scripts_to_execute)
        {
            String command_edited = "python3 "+ fileDir.getAbsolutePath();
            logger.logLine("[DISTRIBUTION-ANALYSIS-PLOTS] Run python script: "+ command_edited);

            Process child = Runtime.getRuntime().exec(command_edited);
            int code = child.waitFor();
            switch (code){
                case 0:
                    break;
                case 1:
                    String message = child.getErrorStream().toString();
                    throw new Exception(message);
            }
        }

        logger.logLine("[DISTRIBUTION-ANALYSIS-PLOTS] Finished comparing background distribution to TF distribution and create plots for outstanding TFs");
    }

    /**
     * creates data needed for distribution plots
     */
    public void perform_distribution_analysis() throws IOException {
        logger.logLine("[DISTRIBUTION-ANALYSIS] Calculate distributions for TFs");

        //search input directories and create output directories
        File f_distr_analysis = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_out_distribution);
        File f_distr_analysis_analysed_tfs = new File(f_distr_analysis.getAbsolutePath()+File.separator+options_intern.folder_out_distribution_analyzed_tfs);
        File f_distr_analysis_analysed_tfs_csv = new File(f_distr_analysis_analysed_tfs.getAbsolutePath()+File.separator+options_intern.file_suffix_distribution_analysis_analysed_tfs);

        File f_out_distr_analysis_tf_tg_scores = new File(f_distr_analysis.getAbsolutePath()+File.separator+options_intern.folder_out_distribution_tf_tg_scores);
        f_out_distr_analysis_tf_tg_scores.mkdir();
        File f_out_distr_analysis_tf_tg_scores_background_distr = new File(f_out_distr_analysis_tf_tg_scores.getAbsolutePath()+File.separator+options_intern.folder_out_distribution_tf_tg_scores_background_distr);
        f_out_distr_analysis_tf_tg_scores_background_distr.mkdir();
        File f_out_distr_analysis_tf_tg_scores_background_distr_ALL = new File(f_out_distr_analysis_tf_tg_scores_background_distr.getAbsolutePath()+File.separator+options_intern.folder_out_distribution_tf_tg_scores_background_distr_ALL);
        f_out_distr_analysis_tf_tg_scores_background_distr_ALL.mkdir();
        File f_out_distr_analysis_tf_tg_scores_background_distr_HM = new File(f_out_distr_analysis_tf_tg_scores_background_distr.getAbsolutePath()+File.separator+options_intern.folder_out_distribution_tf_tg_scores_background_distr_HM);
        f_out_distr_analysis_tf_tg_scores_background_distr_HM.mkdir();
        File f_out_distr_analysis_tf_tg_scores_tf_distr = new File(f_out_distr_analysis_tf_tg_scores.getAbsolutePath()+File.separator+options_intern.folder_out_distribution_tf_tg_scores_tf_distributions);
        f_out_distr_analysis_tf_tg_scores_tf_distr.mkdir();
        File f_out_distr_analysis_tf_tg_scores_tf_distr_ALL = new File(f_out_distr_analysis_tf_tg_scores_tf_distr.getAbsolutePath()+File.separator+options_intern.folder_out_distribution_tf_tg_scores_tf_distributions_ALL);
        f_out_distr_analysis_tf_tg_scores_tf_distr_ALL.mkdir();
        File f_out_distr_analysis_tf_tg_scores_tf_distr_HM = new File(f_out_distr_analysis_tf_tg_scores_tf_distr.getAbsolutePath()+File.separator+options_intern.folder_out_distribution_tf_tg_scores_tf_distributions_HM);
        f_out_distr_analysis_tf_tg_scores_tf_distr_HM.mkdir();

        HashMap<String,HashSet<String>> distinct_hms_tfs = new HashMap<>();
        HashMap<String,HashSet<String>> distinct_tfs_hms = new HashMap<>();

        HashMap<String,String> ensg_symbol = new HashMap<>();
        HashMap<String,HashSet<String>> symbol_ensg = new HashMap<>();

        HashMap<String,HashMap<String,Double>> timepoint_ensg_gene_counts = new HashMap<>();
        HashMap<String,HashMap<String,Double>> group_clash_diff_gene_expression = new HashMap<>();

        HashMap<String,String> composed_tfs = new HashMap<>();
        HashMap<String,HashSet<String>> composed_tfs_tfs = new HashMap<>();

        BufferedReader br_composed_tfs = new BufferedReader(new FileReader(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tepic_postprocessing+File.separator+options_intern.folder_name_tepic_postprocessing_tfs+File.separator+options_intern.file_suffix_tepic_postprocessing_tfs_tfs));
        String line_composed_tfs = "";
        while((line_composed_tfs=br_composed_tfs.readLine())!=null)
        {
            String[] split = line_composed_tfs.split("\t");

            HashSet temp_set = new HashSet();
            for(int i = 1; i < split.length;i++)
            {
                composed_tfs.put(split[i],split[0]);
                temp_set.add(split[i]);
            }
            composed_tfs_tfs.put(split[0],temp_set);
        }
        br_composed_tfs.close();

        BufferedReader br_distinct_hms_tfs = new BufferedReader(new FileReader(f_distr_analysis_analysed_tfs_csv));
        String line_distinct_hms_tfs = br_distinct_hms_tfs.readLine();
        while((line_distinct_hms_tfs=br_distinct_hms_tfs.readLine())!=null)
        {
            String[] split = line_distinct_hms_tfs.split("\t");
            HashSet<String> current_hm_tfs;

            if(distinct_hms_tfs.containsKey(split[1]))
            {
                current_hm_tfs = distinct_hms_tfs.get(split[1]);
            }
            else
            {
                current_hm_tfs=new HashSet<>();

                File f_out_hms_background = new File(f_out_distr_analysis_tf_tg_scores_background_distr_HM.getAbsolutePath()+File.separator+split[1]);
                f_out_hms_background.mkdir();

                File f_out_hms_tf = new File(f_out_distr_analysis_tf_tg_scores_tf_distr_HM.getAbsolutePath()+File.separator+split[1]);
                f_out_hms_tf.mkdir();

            }

            current_hm_tfs.add(split[0]);
            distinct_hms_tfs.put(split[1],current_hm_tfs);

            HashSet<String> current_tf_hms;
            if(distinct_tfs_hms.containsKey(split[0]))
            {
                current_tf_hms = distinct_tfs_hms.get(split[0]);
            }
            else
            {
                current_tf_hms = new HashSet<>();
            }

            current_tf_hms.add(split[1]);
            distinct_tfs_hms.put(split[0],current_tf_hms);
        }
        br_distinct_hms_tfs.close();

        BufferedReader br_ensg_symbol = new BufferedReader(new FileReader(new File(options_intern.tepic_ensg_symbol)));
        String line_ensg_symbol = br_ensg_symbol.readLine();
        while((line_ensg_symbol = br_ensg_symbol.readLine())!=null)
        {
            String[] split=line_ensg_symbol.split("\t");
            if(split.length>1)
            {
                String symbol = split[1].toUpperCase();

                if(composed_tfs.containsKey(symbol))
                {
                    symbol = composed_tfs.get(symbol);
                }

                ensg_symbol.put(split[0].toUpperCase(),symbol);

                HashSet<String> current_ensgs;
                if(symbol_ensg.containsKey(symbol))
                {
                    current_ensgs=symbol_ensg.get(symbol);
                }
                else
                {
                    current_ensgs=new HashSet<>();
                }
                current_ensgs.add(split[0].toUpperCase());
                symbol_ensg.put(symbol,current_ensgs);

            }
        }
        br_ensg_symbol.close();

        File f_genecounts_root = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_deseq2_preprocessing+File.separator+options_intern.folder_name_deseq2_preprocessing_gene_symbols);
        for(File fileDir:f_genecounts_root.listFiles())
        {
            if(fileDir.isFile())
            {
                String name_tp = fileDir.getName().split("\\.")[0];

                HashMap<String,Double> ensg_genecount = new HashMap<>();

                BufferedReader br_gene_counts = new BufferedReader(new FileReader(fileDir));
                String line_gene_counts = br_gene_counts.readLine();
                while((line_gene_counts=br_gene_counts.readLine())!=null)
                {
                    String[] split = line_gene_counts.split("\t");
                    ensg_genecount.put(split[1].toUpperCase(),Double.parseDouble(split[2]));
                }
                br_gene_counts.close();

                timepoint_ensg_gene_counts.put(name_tp,ensg_genecount);
            }
        }

        File f_diff_gene_expr_root = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_deseq2_output);
        for(File fileDir:f_diff_gene_expr_root.listFiles())
        {
            if(fileDir.isFile())
            {
                String[] split_name =fileDir.getName().split("\\.")[0].split("_");
                String name_group_clash = split_name[0]+"_"+split_name[1];

                HashMap<String,Double> current_diff_expr = new HashMap<>();

                BufferedReader br_diff_gene_expr = new BufferedReader(new FileReader(fileDir));
                String line_diff_gene_expr = br_diff_gene_expr.readLine();
                while((line_diff_gene_expr=br_diff_gene_expr.readLine())!=null)
                {
                    String[] split = line_diff_gene_expr.split("\t");
                    current_diff_expr.put(split[0],Double.parseDouble(split[1]));
                }
                br_diff_gene_expr.close();

                group_clash_diff_gene_expression.put(name_group_clash,current_diff_expr);
            }
        }

        for(String key_tf : distinct_tfs_hms.keySet())
        {
            String name_anhaengsel = "";
            if(symbol_ensg.containsKey(key_tf))
            {
                HashSet<String> ensgs = symbol_ensg.get(key_tf);

                HashMap<String,Double> tp_gene_count = new HashMap<>();
                HashMap<String,Double> group_clash_diff_gene_expr = new HashMap<>();

                for(String k_tp : timepoint_ensg_gene_counts.keySet())
                {
                    int gene_count_total = 0;
                    int found_gene_counts = 0;

                    HashMap<String,Double> tp_specific_genecounts = timepoint_ensg_gene_counts.get(k_tp);

                    for(String k_ensgs: ensgs)
                    {
                        if(tp_specific_genecounts.containsKey(k_ensgs))
                        {
                            gene_count_total+=tp_specific_genecounts.get(k_ensgs);
                            found_gene_counts++;
                        }
                    }

                    gene_count_total/=found_gene_counts;

                    tp_gene_count.put(k_tp, (double) gene_count_total);
                }

                for(String k_group_clash : group_clash_diff_gene_expression.keySet())
                {

                    double diff_total = 0;
                    int diff_count = 0;

                    HashMap<String,Double> tp_specific_genecounts = group_clash_diff_gene_expression.get(k_group_clash);

                    for(String k_ensgs: ensgs)
                    {
                        if(tp_specific_genecounts.containsKey(k_ensgs))
                        {
                            diff_total+=tp_specific_genecounts.get(k_ensgs);
                            diff_count++;
                        }
                    }

                    diff_total/=diff_count;

                    group_clash_diff_gene_expr.put(k_group_clash,  diff_total);


                }

                for(String k_group_clash : group_clash_diff_gene_expr.keySet())
                {
                    HashSet<String> available_hms = distinct_tfs_hms.get(key_tf);

                    for (String k_hm : available_hms)
                    {
                        HashSet<String> available_ensgs = new HashSet<>();

                        boolean found_tgene_file = true;
                        if(!options_intern.path_tgen.equals(""))
                        {
                            File f_tgene_input = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tgen+File.separator+options_intern.folder_name_tgen_filter_target_genes+File.separator+k_hm+File.separator+k_group_clash);

                            if(!f_tgene_input.exists())
                            {
                                found_tgene_file=false;
                            }
                            if(found_tgene_file)
                            {
                                for(File fileDir: f_tgene_input.listFiles())
                                {
                                    if(fileDir.isFile())
                                    {
                                        f_tgene_input=fileDir;
                                    }
                                }

                                BufferedReader br_tgene_input = new BufferedReader(new FileReader(f_tgene_input));
                                String line_tgene_input = br_tgene_input.readLine();
                                while((line_tgene_input=br_tgene_input.readLine())!=null)
                                {
                                    String[] split = line_tgene_input.split("\t");
                                    available_ensgs.add(split[0]);
                                }
                                br_tgene_input.close();
                            }

                        }

                        double tf_regression_coefficient = 0.0;

                        File f_hm_tf_coeff_input = new File(options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_put_DYNAMITE + File.separator + k_hm + File.separator+k_group_clash+ File.separator + options_intern.file_suffix_dynamite_output_to_be_plotted);
                        if (f_hm_tf_coeff_input.exists() && f_hm_tf_coeff_input.isFile())
                        {
                            BufferedReader br_hm_tf_coeff = new BufferedReader(new FileReader(f_hm_tf_coeff_input));
                            String line_hm_tf_coeff = br_hm_tf_coeff.readLine();
                            while ((line_hm_tf_coeff = br_hm_tf_coeff.readLine()) != null)
                            {
                                String[] split = line_hm_tf_coeff.split("\t");

                                String[] split_tf_long_name = split[0].split("_");

                                String tf_in_file = split_tf_long_name[0];

                                if(split_tf_long_name.length>1)
                                {
                                    name_anhaengsel=split_tf_long_name[1];
                                }

                                if (tf_in_file.toUpperCase().equals(key_tf.toUpperCase()))
                                {
                                    if (split[1].equals("0"))
                                    {
                                        break;
                                    }
                                    tf_regression_coefficient = Double.parseDouble(split[1]);
                                    break;
                                }

                            }
                        } else {
                            continue;
                        }

                        if(tf_regression_coefficient==0)
                        {
                            continue;
                        }

                        String[] split_clash = k_group_clash.split("_");

                        double current_diff_gene_expr = group_clash_diff_gene_expr.get(k_group_clash);
                        double current_gene_counts = tp_gene_count.get(split_clash[0]) + tp_gene_count.get(split_clash[1]);

                        double current_tf_score = current_gene_counts * current_diff_gene_expr;
                        if (current_tf_score < 0) {
                            current_tf_score *= -1;
                        }

                        File f_group1_hm_tf_target_genes_root = new File(options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_tepic_postprocessing + File.separator + options_intern.folder_name_tepic_postprocessing_output + File.separator + k_hm + File.separator + k_group_clash + File.separator + split_clash[0]);

                        File f_group1_hm_tf_target_genes = new File("");
                        for(File fileDir:f_group1_hm_tf_target_genes_root.listFiles())
                        {
                            if(fileDir.isFile())
                            {
                                f_group1_hm_tf_target_genes=fileDir;
                                break;
                            }
                        }


                        File f_group2_hm_tf_target_genes_root = new File(options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_tepic_postprocessing + File.separator + options_intern.folder_name_tepic_postprocessing_output + File.separator + k_hm + File.separator + k_group_clash + File.separator + split_clash[1]);
                        File f_group2_hm_tf_target_genes = new File("");
                        for(File fileDir:f_group2_hm_tf_target_genes_root.listFiles())
                        {
                            if(fileDir.isFile())
                            {
                                f_group2_hm_tf_target_genes=fileDir;
                                break;
                            }
                        }

                        if (f_group1_hm_tf_target_genes.exists() && f_group1_hm_tf_target_genes.isFile() && f_group2_hm_tf_target_genes.exists() && f_group2_hm_tf_target_genes.isFile())
                        {

                            Map<String,Double> target_genes_scores_group1 = new HashMap<>();
                            Map<String,Double> target_genes_scores_group2 = new HashMap<>();

                            BufferedReader br_group1_hm_tf_target_genes = new BufferedReader(new FileReader(f_group1_hm_tf_target_genes));
                            String line_group1_hm_tf_target_genes=br_group1_hm_tf_target_genes.readLine();

                            String temp_key_tf = key_tf.replace(".",":");
                            String temp_key_tf_2 = key_tf.replace(".",":");
                            if(!name_anhaengsel.equals(""))
                            {
                                temp_key_tf+="_"+name_anhaengsel;
                            }

                            int index_tf_group1_hm_tf_target_genes = -1;
                            String[] split_header_group1_hm_tf_target_gens = line_group1_hm_tf_target_genes.split("\t");
                            for(int i = 0; i < split_header_group1_hm_tf_target_gens.length;i++)
                            {
                                if(split_header_group1_hm_tf_target_gens[i].toUpperCase().equals(temp_key_tf)||split_header_group1_hm_tf_target_gens[i].toUpperCase().equals(temp_key_tf_2))
                                {
                                    index_tf_group1_hm_tf_target_genes=i;
                                    break;
                                }
                            }

                            while((line_group1_hm_tf_target_genes=br_group1_hm_tf_target_genes.readLine())!=null)
                            {
                                String[] split = line_group1_hm_tf_target_genes.split("\t");
                                target_genes_scores_group1.put(split[0],Double.parseDouble(split[index_tf_group1_hm_tf_target_genes]));
                            }
                            br_group1_hm_tf_target_genes.close();

                            BufferedReader br_group2_hm_tf_target_genes = new BufferedReader(new FileReader(f_group2_hm_tf_target_genes));
                            String line_group2_hm_tf_target_genes=br_group2_hm_tf_target_genes.readLine();

                            int index_tf_group2_hm_tf_target_genes = -1;
                            String[] split_header_group2_hm_tf_target_gens = line_group2_hm_tf_target_genes.split("\t");
                            for(int i = 0; i < split_header_group2_hm_tf_target_gens.length;i++)
                            {
                                if(split_header_group2_hm_tf_target_gens[i].toUpperCase().equals(temp_key_tf)||split_header_group2_hm_tf_target_gens[i].toUpperCase().equals(temp_key_tf_2))
                                {
                                    index_tf_group2_hm_tf_target_genes=i;
                                    break;
                                }
                            }

                            while((line_group2_hm_tf_target_genes=br_group2_hm_tf_target_genes.readLine())!=null)
                            {
                                String[] split = line_group2_hm_tf_target_genes.split("\t");
                                target_genes_scores_group2.put(split[0],Double.parseDouble(split[index_tf_group2_hm_tf_target_genes]));
                            }
                            br_group2_hm_tf_target_genes.close();


                            BufferedWriter bw_background_ALL;
                            File f_out_background_ALL = new File(f_out_distr_analysis_tf_tg_scores_background_distr_ALL.getAbsolutePath()+File.separator+options_intern.file_suffix_distribution_analysis_distributions);
                            if(f_out_background_ALL.exists())
                            {
                                bw_background_ALL = new BufferedWriter(new FileWriter(f_out_background_ALL,true));
                            }
                            else
                            {
                                bw_background_ALL = new BufferedWriter(new FileWriter(f_out_background_ALL));
                                bw_background_ALL.write("##HM\tALL\n");
                                bw_background_ALL.write("##TF\tALL\n");
                                bw_background_ALL.write("TARGET_GENE\tTF_TG_SCORE\tHM\tGROUPS\tTF\tTF_COEFF\n");

                            }

                            File f_out_background_HM = new File(f_out_distr_analysis_tf_tg_scores_background_distr_HM.getAbsolutePath()+File.separator+k_hm+File.separator+options_intern.file_suffix_distribution_analysis_distributions);
                            BufferedWriter bw_background_HM;
                            if(f_out_background_HM.exists())
                            {
                                bw_background_HM = new BufferedWriter(new FileWriter(f_out_background_HM,true));
                            }
                            else
                            {
                                bw_background_HM = new BufferedWriter(new FileWriter(f_out_background_HM));
                                bw_background_HM.write("##HM\t"+k_hm+"\n");
                                bw_background_HM.write("##TF\tALL\n");
                                bw_background_HM.write("TARGET_GENE\tTF_TG_SCORE\tHM\tGROUPS\tTF\tTF_COEFF\n");
                            }

                            File f_out_tf_ALL = new File(f_out_distr_analysis_tf_tg_scores_tf_distr_ALL.getAbsolutePath()+File.separator+key_tf+"_"+options_intern.file_suffix_distribution_analysis_distributions);
                            BufferedWriter bw_tf_ALL;
                            if(f_out_tf_ALL.exists())
                            {
                                bw_tf_ALL = new BufferedWriter(new FileWriter(f_out_tf_ALL,true));
                            }
                            else
                            {
                                bw_tf_ALL = new BufferedWriter(new FileWriter(f_out_tf_ALL));
                                bw_tf_ALL.write("##HM\tALL\n");
                                bw_tf_ALL.write("##TF\t"+key_tf+"\n");
                                bw_tf_ALL.write("TARGET_GENE\tTF_TG_SCORE\tHM\tGROUPS\tTF\tTF_COEFF\n");
                            }

                            BufferedWriter bw_tf_HM;
                            File f_out_tf_HM = new File(f_out_distr_analysis_tf_tg_scores_tf_distr_HM.getAbsolutePath()+File.separator+k_hm+File.separator+key_tf+"_"+options_intern.file_suffix_distribution_analysis_distributions);
                            if(f_out_tf_HM.exists())
                            {
                                bw_tf_HM = new BufferedWriter(new FileWriter(f_out_tf_HM,true));
                            }
                            else
                            {
                                bw_tf_HM = new BufferedWriter(new FileWriter(f_out_tf_HM));
                                bw_tf_HM.write("##HM\t"+k_hm+"\n");
                                bw_tf_HM.write("##TF\t"+key_tf+"\n");
                                bw_tf_HM.write("TARGET_GENE\tTF_TG_SCORE\tHM\tGROUPS\tTF\tTF_COEFF\n");
                            }

                            for(String k_target_gene : target_genes_scores_group1.keySet())
                            {
                                if(!options_intern.path_tgen.equals("")&&found_tgene_file)
                                {
                                    if(!available_ensgs.contains(k_target_gene))
                                    {
                                        continue;
                                    }
                                }

                                if(target_genes_scores_group2.containsKey(k_target_gene))
                                {
                                    double gene_score1 = target_genes_scores_group1.get(k_target_gene);
                                    double gene_score2 = target_genes_scores_group2.get(k_target_gene);

                                    double genecount_1 = 0;
                                    double genecount_2 = 0;

                                    double ensg_diff_gene_expr = 0;
                                    if(timepoint_ensg_gene_counts.get(split_clash[0]).containsKey(k_target_gene))
                                    {
                                        genecount_1=timepoint_ensg_gene_counts.get(split_clash[0]).get(k_target_gene);
                                    }

                                    if(timepoint_ensg_gene_counts.get(split_clash[1]).containsKey(k_target_gene))
                                    {
                                        genecount_2=timepoint_ensg_gene_counts.get(split_clash[1]).get(k_target_gene);
                                    }

                                    if(group_clash_diff_gene_expr.containsKey(k_group_clash))
                                    {
                                        ensg_diff_gene_expr=group_clash_diff_gene_expr.get(k_group_clash);
                                    }

                                    double tg_score_1 = genecount_1*ensg_diff_gene_expr*gene_score1;
                                    double tg_score_2 = genecount_2*ensg_diff_gene_expr*gene_score2;
                                    if(tg_score_1<0)
                                    {
                                        tg_score_1*=-1;
                                    }
                                    if(tg_score_2<0)
                                    {
                                        tg_score_2*=-1;
                                    }

                                    double tg_score_cumm = tg_score_1 + tg_score_2;

                                    double tf_tg_score = (current_tf_score*tg_score_cumm)*tf_regression_coefficient;

                                    if(tf_tg_score<0)
                                    {
                                        tf_tg_score*=-1;
                                    }

                                    StringBuilder sb_line = new StringBuilder();
                                    sb_line.append(k_target_gene);
                                    sb_line.append("\t");
                                    sb_line.append(tf_tg_score);
                                    sb_line.append("\t");
                                    sb_line.append(k_hm);
                                    sb_line.append("\t");
                                    sb_line.append(k_group_clash);
                                    sb_line.append("\t");
                                    sb_line.append(key_tf);
                                    sb_line.append("\t");
                                    sb_line.append(tf_regression_coefficient);

                                    bw_background_HM.write(sb_line.toString());
                                    bw_background_HM.newLine();

                                    bw_background_ALL.write(sb_line.toString());
                                    bw_background_ALL.newLine();

                                    bw_tf_ALL.write(sb_line.toString());
                                    bw_tf_ALL.newLine();

                                    bw_tf_HM.write(sb_line.toString());
                                    bw_tf_HM.newLine();
                                }
                            }
                            bw_background_ALL.close();
                            bw_background_HM.close();
                            bw_tf_ALL.close();
                            bw_tf_HM.close();
                        }
                    }
                }
            }
            else
            {
                continue;
            }
        }


        logger.logLine("[DISTRIBUTION-ANALYSIS] Finished calculating distributions");
    }

    /**
     * creates the HTML report
     * it creates also one important file (available tfs in HMs and which stages) for distribution analysis
     */
    public void create_overview_html_report_before_distribution_analysis() throws Exception {
        logger.logLine("[WEBSITE] Start creating overview website.");

        File f_website_css = new File(options_intern.path_to_COM2POSE+File.separator+"ext"+File.separator+"WEBSITE");

        File f_move_css = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_out_website+File.separator+options_intern.folder_out_website_basics);
        f_move_css.mkdir();

        File f_output_website = new File(options_intern.com2pose_working_directory+ File.separator+options_intern.folder_out_website);
        f_output_website.mkdir();

        File f_output_website_htmls = new File(f_output_website.getAbsolutePath()+File.separator+options_intern.folder_out_website_htmls_regression_coefficients);
        f_output_website_htmls.mkdir();

        String command = "cp -u -r " + f_website_css.getAbsolutePath() + " " + f_move_css.getAbsolutePath();
        Process child = Runtime.getRuntime().exec(command);
        logger.logLine("[WEBSITE] Copy CSS files: " + command);
        int code = child.waitFor();
        switch (code) {
            case 0:
                break;
            case 1:
                String message = child.getErrorStream().toString();
                throw new Exception(message);
        }

        for(int i_dont_know = 0; i_dont_know < 2; i_dont_know++)
        {

            HashMap<String,HashMap<String,HashSet<String>>> distinct_tf_hm_diff_same = new HashMap<>();

            String html_tail = "</body>\n" +
                    "</html>";

            BufferedWriter bw_home = new BufferedWriter(new FileWriter(f_output_website+File.separator+options_intern.html_report_home_regression_coefficient_analysis));

            //bw_home.write(sb_home_front.toString());
            bw_home.write(get_header_html("HOME", options_intern.analysis_types_regression_coefficient_analysis));

            threshold_folders_filled =true;

            for(Double d: options_intern.plot_th_coefficient)
            {
                bw_home.append(write_table_html(d,"HOME"));
            }

            bw_home.write(html_tail);

            bw_home.close();


            StringBuilder sb_parameter = new StringBuilder(get_header_html("HOME",options_intern.analysis_types_regression_coefficient_analysis));

            sb_parameter.append(" <script>\n" +
                    " document.title = \"PARAMETERS\";\n" +
                    " </script>\n");

            /*
            sb_parameter.append("<button class=\"button_expandable\" id=\"BUTTON_NAME\" aria-expanded=\"false\" ondblclick=\"expand_collapse('BUTTON_NAME','TABLE_NAME')\"> FILTER_NAME\n");
            sb_parameter.append("<div style=\"display: none;background-color: white;color:black;\" id=\"TABLE_NAME\">\n");

            sb_parameter.append("\t\t<table style=\"width:100%;font-size:15px;text-align: left;\">\n");

            sb_parameter.append("\t\t\t<tr>");
            sb_parameter.append("\t\t\t\t<th>Parameter</th>");
            sb_parameter.append("\t\t\t\t<th>Value</th>");
            sb_parameter.append("\t\t\t\t<th>Explanation</th>");
            sb_parameter.append("\t\t\t</tr>");

            sb_parameter.append("\t\t</table>");
            sb_parameter.append("</div>\n</button>\n");
             */


            //mix_option_parameters
            {
                sb_parameter.append("<button class=\"button_expandable\" id=\"button_preprocessing_mix_options\" aria-expanded=\"false\" ondblclick=\"expand_collapse('button_preprocessing_mix_options','table_preprocessing_mix_options')\"> Preprocessing Mix Options\n");
                sb_parameter.append("<div style=\"display: none;background-color: white;color:black;\" id=\"table_preprocessing_mix_options\">\n");

                sb_parameter.append("\t\t<table style=\"width:100%;font-size:15px;text-align: left;\">\n");
                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>Parameter</th>");
                sb_parameter.append("\t\t\t\t<th>Value</th>");
                sb_parameter.append("\t\t\t\t<th>Explanation</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>mix_level</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.mix_level + "</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: if set a mix of either the Histone Modification Level (HM_LEVEL) or the Sample Level (SAMPLE_LEVEL) will be performed, mix_option is required\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>mix_option</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.mix_option + "</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: set histone marks or samples will be mixed with option: UNION (all peaks of all HMs will be used), INTERSECTION (only peaks, which are in all HMs will be used)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                if(options_intern.mix_option.equals("INTERSECTION"))
                {
                    sb_parameter.append("\t\t\t<tr>");
                    sb_parameter.append("\t\t\t\t<th>mix_occurence_intersection</th>");
                    sb_parameter.append("\t\t\t\t<th>"+options_intern.mix_occurence_intersection+"</th>");
                    sb_parameter.append("\t\t\t\t<th>#[OPT]: minimal occurence of peaks in intersection (only applied if mix_option is set to INTERSECTION), if set to zero it means it must be in all samples /histone modifications available, default 2\n</th>");
                    sb_parameter.append("\t\t\t</tr>");
                }
                sb_parameter.append("\t\t</table>");
                sb_parameter.append("</div>\n</button>\n");
            }

            //blacklist parameters
            {
                sb_parameter.append("<button class=\"button_expandable\" id=\"button_blacklist_parameters\" aria-expanded=\"false\" ondblclick=\"expand_collapse('button_blacklist_parameters','table_blacklist_parameters')\"> Blacklist Parameters\n");
                sb_parameter.append("<div style=\"display: none;background-color: white;color:black;\" id=\"table_blacklist_parameters\">\n");
                sb_parameter.append("\t\t<table style=\"width:100%;font-size:15px;text-align: left;\">\n");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>Parameter</th>");
                sb_parameter.append("\t\t\t\t<th>Value</th>");
                sb_parameter.append("\t\t\t\t<th>Explanation</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>black_list_dir</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.black_list_dir+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: if blacklist filter should be used provide path to BED file here (BED files can be found at https://github.com/Boyle-Lab/Blacklist/tree/master/lists)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>black_list_signals</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.black_list_signals.toString()+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: which regions should be ignored: Low_Mappability;High_Signal_Region\n</th>");
                sb_parameter.append("\t\t\t</tr>");


                sb_parameter.append("\t\t</table>");
                sb_parameter.append("</div>\n</button>\n");
            }

            //DESeq2 parameters
            {
                sb_parameter.append("<button class=\"button_expandable\" id=\"button_deseq2_parameters\" aria-expanded=\"false\" ondblclick=\"expand_collapse('button_deseq2_parameters','table_deseq2_parameters')\"> DESeq2 Parameters\n");
                sb_parameter.append("<div style=\"display: none;background-color: white;color:black;\" id=\"table_deseq2_parameters\">\n");

                sb_parameter.append("\t\t<table style=\"width:100%;font-size:15px;text-align: left;\">\n");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>Parameter</th>");
                sb_parameter.append("\t\t\t\t<th>Value</th>");
                sb_parameter.append("\t\t\t\t<th>Explanation</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>deseq2_input_directory</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.deseq2_input_directory+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[REQ]: count files in nfcore RNA-seq format (each line is a count of one gene), directory must be ordered like: TP1 - samples_1,...,samples_n;....\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>deseq2_input_gene_id</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.deseq2_input_gene_id+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[REQ]: gene ID file from nfcore RNA-seq (each line names one gene ID (ENSG) - must be same order as deseq2_input_directory files\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>deseq2_biomart_dataset_species</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.deseq2_biomart_dataset_species+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[REQ]: biomart dataset name for the species which is used in RNA-seq and CHIP-seq data, (https://www.bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>deseq2_count_threshold</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.deseq2_count_threshold+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: minimum count over all samples of two timepoints for DESeq2, default: 0\n</th>");
                sb_parameter.append("\t\t\t</tr>");
                sb_parameter.append("\t\t</table>");
                sb_parameter.append("</div>\n</button>\n");
            }

            //TEPIC parameters
            {
                sb_parameter.append("<button class=\"button_expandable\" id=\"button_tepic_parameters\" aria-expanded=\"false\" ondblclick=\"expand_collapse('button_tepic_parameters','table_tepic_parameters')\"> TEPIC Parameters\n");
                sb_parameter.append("<div style=\"display: none;background-color: white;color:black;\" id=\"table_tepic_parameters\">\n");

                sb_parameter.append("\t\t<table style=\"width:100%;font-size:15px;text-align: left;\">\n");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>Parameter</th>");
                sb_parameter.append("\t\t\t\t<th>Value</th>");
                sb_parameter.append("\t\t\t\t<th>Explanation</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_input_directory [ORIGINAL]</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tepic_input_original+"</th>");
                sb_parameter.append("\t\t\t\t<th></th>");
                sb_parameter.append("\t\t\t</tr>");

                if(!options_intern.tepic_input_prev.equals(""))
                {
                    sb_parameter.append("\t\t\t<tr>");
                    sb_parameter.append("\t\t\t\t<th>tepic_input_directory [LAST FILTER]</th>");
                    sb_parameter.append("\t\t\t\t<th>"+options_intern.tepic_input_directory+"</th>");
                    sb_parameter.append("\t\t\t\t<th></th>");
                    sb_parameter.append("\t\t\t</tr>");
                }

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_input_ref_genome</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tepic_input_ref_genome+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[REQ]:fasta files (in RefSeq format, without \\\"chr\\\" prefix)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_path_pwms</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tepic_path_pwms+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[REQ]: path to position specific energy matrix used for TRAP (different matrices can be found in ~/COM2POSE/ext/TEPIC/TEPIC/PWMs/2.1)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_cores</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tepic_cores+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: number of cores\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_bed_chr_sign</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tepic_bed_chr_sign+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: bedgraph file containing open chromatin signal, e.g. DNase1-seq, or Histone-Mark signal\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_column_bedfile</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tepic_column_bedfile+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: column in the tepic_input_directory file containing the average per base signal within a peak. If this option is used, the tepic_bed_chr_sign option must not be used\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_gene_annot</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tepic_gene_annot+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[REQ]: gene annotation file, required to generate the gene view, required for TPM filter\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_window_size</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tepic_window_size+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: size of the window to be considered to generate gene view (default 50kb)]\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_onlyDNasePeaks</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tepic_onlyDNasePeaks+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: path to annotation file and annotate only DNase peaks that are within a window specified by the tepic_window_size option around all genes contained in the gene annotation file specified by this option\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_exponential_decay</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tepic_exponential_decay+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: indicating whether exponential decay should be used (default TRUE)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_not_norm_peak_length</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tepic_not_norm_peak_length+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: flag to be set if affinities should not be normalised by peak length\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_not_generated</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tepic_not_generated+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: flag to be set if peak features for peak length and peak counts should not be generated\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_original_decay</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tepic_original_decay+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: if tepic_bed_chr_sign or tepic_column_bedfile is used together with this flag, the original (Decay-)Scaling formulation of TEPIC is used to compute gene-TF scores\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_psems_length</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tepic_psems_length+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: path to a tab delimited file containing the length of the used PSEMs\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_entire_gene_body</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tepic_entire_gene_body+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: flag to be set if the entire gene body should be screened for TF binding. The search window is extended by a region half of the size that is specified by the tepic_window_size option upstream of the genes 5' TSS\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_zipped</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tepic_zipped+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: flag indicating that the output of TEPIC should be zipped\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_background_seq</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tepic_background_seq+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: path to a set of background sequences that should be used to compute to generate a binary score for TF binding. Mutually exclusive to the tepic_2bit option\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_2bit</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tepic_2bit+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: path to a 2bit representation of the reference genome, required to generate a binary score for TF binding. The binary score is generated in addition to the standard affinity values. Mutually exclusive to the tepic_background_seq option\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_pvalue</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tepic_pvalue+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: p-value cut off used to determine a cut off to derive a binary score for TF binding (default 0.05)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_minutes_per_chr</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tepic_minutes_per_chr+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: minutes that should be spend at most per chromosome to find matching random regions (default 3)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_chr_prefix</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tepic_chr_prefix+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: flag indicating that the reference genome contains a chr prefix\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_transcript_based</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tepic_transcript_based+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: flag indicating that the annotation should be transcript and not gene based\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_loop_list</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tepic_loop_list+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: loop list file containing chromatin contacts\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_loop_windows</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tepic_loop_windows+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: size of the loop window used around a genes promoter to link chromatin loops to genes (default 5000)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_only_peak_features</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tepic_only_peak_features+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: parameter to be set if only peak features should be computed (default FALSE)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_tpm_cutoff</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tepic_tpm_cutoff+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: set (T)ranscripts (P)er (M)illion cutoff, default: TPM filter not active\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_ensg_symbol</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tepic_ensg_symbol+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[REQ]: IF NOT SET MAPPING WILL BE PERFORMED AUTOMATICALLY - path to input file ensg to gene symbol file\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t</table>");
                sb_parameter.append("</div>\n</button>\n");

            }

            //TGEN parameters
            if(!options_intern.path_tgen.equals(""))
            {
                sb_parameter.append("<button class=\"button_expandable\" id=\"button_tgen_parameters\" aria-expanded=\"false\" ondblclick=\"expand_collapse('button_tgen_parameters','table_tgen_parameters')\"> TGene Parameters\n");
                sb_parameter.append("<div style=\"display: none;background-color: white;color:black;\" id=\"table_tgen_parameters\">\n");

                sb_parameter.append("\t\t<table style=\"width:100%;font-size:15px;text-align: left;\">\n");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>Parameter</th>");
                sb_parameter.append("\t\t\t\t<th>Value</th>");
                sb_parameter.append("\t\t\t\t<th>Explanation</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tgen_no_closest_locus</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tgen_no_closest_locus+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: if no locus is found within window size, the nearest locus is used, default:false, meaning locus is used\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tgen_no_closest_tss</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tgen_no_closest_tss+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: if no tss is found within window size, the nearest locus is used, default:false, meaning locus is used\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tgen_max_link_distances</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tgen_max_link_distances+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: window size of tgen, default: 500000\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tgen_pvalue</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tgen_pvalue+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: max pvalue, which is accepted, default 0.05\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tgen_self_regulatory</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.tgen_self_regulatory+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[REQ]: should self regulatory TFs {OPT-SELF-REG} be increased? default: false\n</th>");
                sb_parameter.append("\t\t\t</tr>");
                if(options_intern.tgen_self_regulatory)
                {
                    sb_parameter.append("\t\t\t<tr>");
                    sb_parameter.append("\t\t\t\t<th>tgen_consensus</th>");
                    sb_parameter.append("\t\t\t\t<th>"+options_intern.tgen_consensus+"</th>");
                    sb_parameter.append("\t\t\t\t<th>#[OPT]: {OPT-SELF-REG} value for consus approach (e.g. 0.5 = 50:50 TGene:TEPIC, 0.4 = 40:60 TGene:TEPIC), default: 0.5\n</th>");
                    sb_parameter.append("\t\t\t</tr>");

                    sb_parameter.append("\t\t\t<tr>");
                    sb_parameter.append("\t\t\t\t<th>tgen_consensus_calc</th>");
                    sb_parameter.append("\t\t\t\t<th>"+options_intern.tgen_consensus_calc+"</th>");
                    sb_parameter.append("\t\t\t\t<th>#[OPT]: {OPT-SELF-REG} can be \"INCREASE_TGENE_TFS\" (increases TEPIC TF affinities for TFs found in TGENE at target gene) or \"DECREASE_NOT_TGENE_TFs (decreases TEPIC TF affinities for TFs not found in TGENE at target gene)\".\n</th>");
                    sb_parameter.append("\t\t\t</tr>");

                    sb_parameter.append("\t\t\t<tr>");
                    sb_parameter.append("\t\t\t\t<th>tgen_mt_writing</th>");
                    sb_parameter.append("\t\t\t\t<th>"+options_intern.tgen_mt_writing+"</th>");
                    sb_parameter.append("\t\t\t\t<th>#[OPT]: {OPT-SELF-REG} if TGENE consensus is used please specify the writing of the Mitochondrial DNA chromosome in Peak Files, default: MT\n</th>");
                    sb_parameter.append("\t\t\t</tr>");

                }
                sb_parameter.append("\t\t</table>");
                sb_parameter.append("</div>\n</button>\n");
            }

            //DYNAMITE parameters
            {
                sb_parameter.append("<button class=\"button_expandable\" id=\"button_dynamite_parameters\" aria-expanded=\"false\" ondblclick=\"expand_collapse('button_dynamite_parameters','table_dynamite_parameters')\"> DYNAMITE parameters\n");
                sb_parameter.append("<div style=\"display: none;background-color: white;color:black;\" id=\"table_dynamite_parameters\">\n");

                sb_parameter.append("\t\t<table style=\"width:100%;font-size:15px;text-align: left;\">\n");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>Parameter</th>");
                sb_parameter.append("\t\t\t\t<th>Value</th>");
                sb_parameter.append("\t\t\t\t<th>Explanation</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>dynamite_preprocessing_integrate_data_consider_geneFile</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.dynamite_preprocessing_integrate_data_consider_geneFile+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: File containing gene IDs that should be considered\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>dynamite_out_var</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.dynamite_out_var+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[REQ]: Name of the response variable default: Expression (in this pipeline)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>dynamite_cores</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.dynamite_cores+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: Number of the cores to use (1 as default)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>dynamite_alpha</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.dynamite_alpha+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: Alpha parameter stepsize (0.1 as default)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>dynamite_testsize</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.dynamite_testsize+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]:Size of test data (0.2 as default)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>dynamite_Ofolds</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.dynamite_Ofolds+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: Number of outer folds for model validation (3 as default)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>dynamite_Ifolds</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.dynamite_Ifolds+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: Number of inner cross validation folds (6 as default)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>dynamite_balanced</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.dynamite_balanced+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: Flag indicating whether the data should be balanced through downsampling (default TRUE)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>dynamite_performance</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.dynamite_performance+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: Flag indicating whether performance measures should be computed (default TRUE)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>dynamite_randomise</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.dynamite_randomise+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: Flag indicating whether a model should be learned on randomised data (default FALSE)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t</table>");
                sb_parameter.append("</div>\n</button>\n");

            }

            //PLOT parameters
            {
                sb_parameter.append("<button class=\"button_expandable\" id=\"button_plot_parameters\" aria-expanded=\"false\" ondblclick=\"expand_collapse('button_plot_parameters','table_plot_parameters')\"> PLOT Parameters\n");
                sb_parameter.append("<div style=\"display: none;background-color: white;color:black;\" id=\"table_plot_parameters\">\n");

                sb_parameter.append("\t\t<table style=\"width:100%;font-size:15px;text-align: left;\">\n");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>Parameter</th>");
                sb_parameter.append("\t\t\t\t<th>Value</th>");
                sb_parameter.append("\t\t\t\t<th>Explanation</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>plot_th_coefficient</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.plot_th_coefficient.toString()+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: thresholds for coefficent plot and overall plots, for each threshold it creates a plot for each timepoint and an overall plot for each HistoneMod\n" +
                        "#e.g. 0.1 means it creates a plot of the coefficient range ([-1;1]), it uses all TFs of [-1,-0.1] and [0.1,1]\n" +
                        "#default: 0.1;0.2;0.3;0.4;0.5;0.6</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>plot_cutoff_tps</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.plot_cutoff_tps+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: in how many TPs should a TF be found, default: 2\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>plot_cutoff_hms</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.plot_cutoff_hms+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: in how many HMs should a TF be found, default: 2\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>plot_cutoff_gcs</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.plot_cutoff_gcs+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: minimum gene counts, default: 100\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>plot_top_k_genes</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.plot_top_k_genes+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: top k target genes for TFs, default: 30\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t</table>");
                sb_parameter.append("</div>\n</button>\n");

            }

            //HTML Report PARAMETERS
            {
                sb_parameter.append("<button class=\"button_expandable\" id=\"button_html_report_parameters\" aria-expanded=\"false\" ondblclick=\"expand_collapse('button_html_report_parameters','table_html_report_parameters')\"> HTML Report Parameters\n");
                sb_parameter.append("<div style=\"display: none;background-color: white;color:black;\" id=\"table_html_report_parameters\">\n");

                sb_parameter.append("\t\t<table style=\"width:100%;font-size:15px;text-align: left;\">\n");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>Parameter</th>");
                sb_parameter.append("\t\t\t\t<th>Value</th>");
                sb_parameter.append("\t\t\t\t<th>Explanation</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>html_report_interesting_tfs</th>");
                sb_parameter.append("\t\t\t\t<th>"+options_intern.website_interesting_tfs+"</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: list of TFs which you are interested in - is only used to search fast for known TFs in the results, it does not affect results\n" +
                        "#e.g. website_interesting_tfs=\"STAT3;GATA3;NFIB\"</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t</table>");
                sb_parameter.append("</div>\n</button>\n");
            }


            /*
            sb_parameter.append("\t\t\t<tr>");
            sb_parameter.append("\t\t\t\t<th></th>");
            sb_parameter.append("\t\t\t\t<th></th>");
            sb_parameter.append("\t\t\t\t<th></th>");
            sb_parameter.append("\t\t\t</tr>");
             */


            sb_parameter.append(html_tail);

            BufferedWriter bw_parameter = new BufferedWriter(new FileWriter(new File(f_output_website+File.separator+"PARAMETERS.html")));
            bw_parameter.write(sb_parameter.toString());
            bw_parameter.close();


            HashMap<String,ArrayList<String>> tf_gene_count = new HashMap<>();

            for(File fileDir: threshold_folders)
            {
                HashSet<String> tfs_to_create_pages = new HashSet<>();
                HashSet<String> possible_hms = new HashSet<>();

                File f_interactive_plots_root = new File(f_output_website.getAbsolutePath()+ File.separator+options_intern.folder_out_website_interactive_plots+File.separator+fileDir.getName());

                StringBuilder sb_threshold = new StringBuilder(get_header_html("THRESHOLD",options_intern.analysis_types_regression_coefficient_analysis));
                sb_threshold.append(" <script>\n" +
                        " document.title = \"TH: "+fileDir.getName()+"\";\n" +
                        " </script>\n");
                sb_threshold.append("<div class=\"w3-row-padding w3-padding-64 w3-container\"><div class=\"w3-content\"><h1>Current threshold: "+fileDir.getName()+"</h1></div>\n</div>\n");


                HashSet<String> total_number_tfs = new HashSet<>();

                for(File fileDir_hm : f_interactive_plots_root.listFiles())
                {
                    HashSet<String> total_number_tfs_hm = new HashSet<>();

                    possible_hms.add(fileDir_hm.getName());

                    sb_threshold.append("<div class=\"w3-row-padding w3-padding-64 w3-container\">\n" +
                            "  <div class=\"w3-content\">\n" +
                            "    <div class=\"w3-twothird\">\n" +
                            "      <h1>");
                    sb_threshold.append(fileDir_hm.getName());
                    sb_threshold.append("</h1>\n" +
                            "\t  \n" +
                            "\t <h4> Different timepoints / conditions </h4> \n" +
                            "\t  \n");
                    sb_threshold.append("</div>\n");
                    sb_threshold.append("</div>\n");
                    sb_threshold.append("<div class=\"container\" style=\"\n" +
                            "    width: 100%;\n" +
                            "    position: relative;\n" +
                            "    display: block;\n" +
                            "    overflow-x: scroll;overflow-y: scroll;\">\n" +
                            "\t  <iframe id=\"igraph_"+fileDir_hm.getName()+"_threshold_"+fileDir.getName()+"_different_stages.html"+"\" scrolling=\"no\" style=\"border:none;\" seamless=\"seamless\" src=\"");
                    String relative_path_diff_th = ".."+File.separator+".."+File.separator+options_intern.folder_out_website_interactive_plots+File.separator+fileDir.getName()+File.separator+fileDir_hm.getName()+File.separator+options_intern.folder_out_website_interactive_plots_overview +File.separator+ fileDir_hm.getName()+"_threshold_"+fileDir.getName()+"_different_stages.html";

                    sb_threshold.append(relative_path_diff_th);
                    sb_threshold.append("\" height=\"400\" width=\"2500\" overflow=\"scroll\"></iframe>\n");
                    sb_threshold.append("</div>\n");

                    sb_threshold.append("    <div class=\"w3-twothird\">\n");
                    sb_threshold.append("  <div class=\"w3-content\">\n");

                    //gene count table for plot

                    //collapse thing
                    sb_threshold.append("<button class=\"button_expandable\" id=\"button_"+fileDir_hm.getName()+"_diff\" aria-expanded=\"false\" ondblclick=\"expand_collapse('button_"+fileDir_hm.getName()+"_diff','table_"+fileDir_hm.getName()+"_diff')\"> GeneCounts\n");
                    //tooltip
                    sb_threshold.append("<div class=\"tooltip\"><img src=\".."+File.separator+".."+File.separator+options_intern.folder_out_website_basics+File.separator+options_intern.folder_out_website_basics_website+File.separator+options_intern.folder_out_website_basics_website_images+File.separator+"information.png"+"\" style=\"width:45px;height:40px;\"/>"+
                            "  <span class=\"tooltiptext\">GeneCount is the unprocessed data counted by nfcore RNA-seq. A TF which is at least expressed with GeneCount 1000 is considered highly expressed and active.</span>\n" +
                            "</div>");
                    sb_threshold.append("<div style=\"display: none;background-color: white;color:black;\" id=\"table_"+fileDir_hm.getName()+"_diff\">\n");

                    HashSet<String> total_numbers_tfs_hm_diff = new HashSet<>();
                    sb_threshold.append("<h4>Gene Count threshold: "+options_intern.plot_cutoff_gcs+"</h4>");
                    sb_threshold.append("<h5><i>Click on TF for detailed information - if no Button is available it means that this TF was eliminated by a filter.</i></h5>\n");
                    sb_threshold.append("\t\t<table style=\"width:100%;font-size:15px;\">\n");

                    File f_gene_counts_input_different = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_out_analysis_data+File.separator+options_intern.folder_out_analysis_data_TP_LEVEL+File.separator+fileDir_hm.getName()+File.separator+fileDir.getName()+File.separator+options_intern.file_suffix_analysis_plot_data_hm_level_different);
                    BufferedReader br_gc_different = new BufferedReader(new FileReader(f_gene_counts_input_different));
                    String line_gc_different ="";
                    while((line_gc_different=br_gc_different.readLine())!=null)
                    {
                        String[] split = line_gc_different.split("\t");
                        sb_threshold.append("\t\t\t<tr>\n");

                        tfs_to_create_pages.add(split[0]);

                        String tf_key ="";
                        ArrayList<String> tf_key_set = new ArrayList<>();

                        int count = 0;
                        for(String s: split)
                        {
                            if(count==0)
                            {
                                tf_key=s;
                            }
                            if(count!=0)
                            {
                                tf_key_set.add(s);
                            }
                            sb_threshold.append("\t\t\t\t<th>");
                            if(count!=0 || s.equals("TF"))
                            {
                                sb_threshold.append(s.toUpperCase());
                            }
                            else
                            {
                                File f_try = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_out_website+File.separator+options_intern.folder_out_website_htmls_regression_coefficients+File.separator+fileDir.getName()+File.separator+options_intern.folder_out_website_htmls_TFs+File.separator+s.toUpperCase()+".html");
                                if(f_try.exists())
                                {
                                    sb_threshold.append("<a href='");
                                    //sb_threshold.append(f_try.getAbsolutePath());
                                    sb_threshold.append("TFs"+File.separator+f_try.getName());
                                    sb_threshold.append("' target='_blank'><button class=\"button\">"+s.toUpperCase()+"</button>");
                                    sb_threshold.append("</a>");

                                    total_number_tfs.add(s.toUpperCase());
                                    total_number_tfs_hm.add(s.toUpperCase());
                                    total_numbers_tfs_hm_diff.add(s.toUpperCase());

                                    HashMap<String,HashSet<String>> current_tf_hm;

                                    if(distinct_tf_hm_diff_same.containsKey(s.toUpperCase()))
                                    {
                                        current_tf_hm=distinct_tf_hm_diff_same.get(s.toUpperCase());
                                    }
                                    else
                                    {
                                        current_tf_hm = new HashMap<>();
                                    }

                                    HashSet<String> current_tf_hm_stage = new HashSet<>();

                                    if(current_tf_hm.containsKey(fileDir_hm.getName()))
                                    {
                                        current_tf_hm_stage = current_tf_hm.get(fileDir_hm.getName());
                                    }
                                    else
                                    {
                                        current_tf_hm_stage= new HashSet<>();
                                    }

                                    current_tf_hm_stage.add("DIFFERENT_TPS");
                                    current_tf_hm.put(fileDir_hm.getName(),current_tf_hm_stage);
                                    distinct_tf_hm_diff_same.put(s.toUpperCase(),current_tf_hm);

                                }
                                else
                                {
                                    sb_threshold.append(s.toUpperCase());
                                }

                            }
                            sb_threshold.append("\t\t\t\t</th>\n");
                            count++;
                        }

                        tf_gene_count.put(tf_key.toUpperCase(),tf_key_set);

                        sb_threshold.append("\t\t\t</tr>\n");

                    }
                    br_gc_different.close();

                    sb_threshold.append("\t\t</table>\n");

                    //collapse thing
                    sb_threshold.append("</div>\n");
                    sb_threshold.append("</button>\n");

                    sb_threshold.append("<h5> A total number of "+total_numbers_tfs_hm_diff.size()+" distinct TFs are considered in different time points group </h5>\n");

                    HashSet<String> total_numbers_tfs_hm_same = new HashSet<>();
                    sb_threshold.append("\t <h4> Same timepoints / conditions </h4> \n" +
                            "\t  \n" );
                    sb_threshold.append("</div>\n");
                    sb_threshold.append("</div>\n");
                    sb_threshold.append("<div class=\"container\" style=\"\n" +
                            "    width: 100%;\n" +
                            "    position: relative;\n" +
                            "    display: block;\n" +
                            "    overflow-x: scroll;overflow-y: scroll;\">\n" +
                            "\t  <iframe id=\"igraph_"+fileDir_hm.getName()+"_threshold_"+fileDir.getName()+"_same_stages.html"+"\" scrolling=\"no\" style=\"border:none;\" seamless=\"seamless\" src=\"");
                    String relative_path_same_th = ".."+File.separator+".."+File.separator+options_intern.folder_out_website_interactive_plots+File.separator+fileDir.getName()+File.separator+fileDir_hm.getName()+File.separator+options_intern.folder_out_website_interactive_plots_overview +File.separator+ fileDir_hm.getName()+"_threshold_"+fileDir.getName()+"_same_stages.html";
                    sb_threshold.append(relative_path_same_th);
                    sb_threshold.append("\" height=\"400\" width=\"2500\"></iframe>\n");
                    sb_threshold.append("</div>\n");

                    sb_threshold.append("    <div class=\"w3-twothird\">\n");
                    sb_threshold.append("  <div class=\"w3-content\">\n");


                    //collapse thing
                    sb_threshold.append("<button class=\"button_expandable\" id=\"button_"+fileDir_hm.getName()+"_same\" aria-expanded=\"false\" ondblclick=\"expand_collapse('button_"+fileDir_hm.getName()+"_diff','table_"+fileDir_hm.getName()+"_same')\"> GeneCounts\n");
                    //tooltip
                    sb_threshold.append("<div class=\"tooltip\"><img src=\".."+File.separator+".."+File.separator+options_intern.folder_out_website_basics+File.separator+options_intern.folder_out_website_basics_website+File.separator+options_intern.folder_out_website_basics_website_images+File.separator+"information.png"+"\" style=\"width:45px;height:40px;\"/>"+
                            "  <span class=\"tooltiptext\">GeneCount is the unprocessed data counted by nfcore RNA-seq. A TF which is at least expressed with GeneCount 1000 is considered highly expressed and active.</span>\n" +
                            "</div>");
                    sb_threshold.append("<div style=\"display: none;background-color: white;color:black;\" id=\"table_"+fileDir_hm.getName()+"_same\">\n");

                    //gene count table for plot
                    sb_threshold.append("<h4>Gene Count threshold: "+options_intern.plot_cutoff_gcs+"</h4>");
                    sb_threshold.append("<h5><i>Click on TF for detailed information - if no Button is available it means that this TF was eliminated by a filter.</i></h5>\n");
                    sb_threshold.append("\t\t<table style=\"width:100%;font-size:15px;\">\n");

                    File f_gene_counts_input_same = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_out_analysis_data+File.separator+options_intern.folder_out_analysis_data_TP_LEVEL+File.separator+fileDir_hm.getName()+File.separator+fileDir.getName()+File.separator+options_intern.file_suffix_analysis_plot_data_hm_level_same);
                    BufferedReader br_gc_same = new BufferedReader(new FileReader(f_gene_counts_input_same));
                    String line_gc_same ="";
                    while((line_gc_same=br_gc_same.readLine())!=null)
                    {
                        String[] split = line_gc_same.split("\t");
                        sb_threshold.append("\t\t\t<tr>\n");

                        tfs_to_create_pages.add(split[0]);

                        String tf_key="";
                        ArrayList<String> tf_key_set = new ArrayList<>();

                        int count = 0;
                        for(String s: split)
                        {
                            if(count==0)
                            {
                                tf_key=s;
                            }
                            if(count!=0)
                            {
                                tf_key_set.add(s);
                            }

                            sb_threshold.append("\t\t\t\t<th>");
                            if(count!=0 || s.equals("TF"))
                            {
                                sb_threshold.append(s.toUpperCase());
                            }
                            else
                            {
                                File f_try = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_out_website+File.separator+options_intern.folder_out_website_htmls_regression_coefficients+File.separator+fileDir.getName()+File.separator+options_intern.folder_out_website_htmls_TFs+File.separator+s.toUpperCase()+".html");
                                if(f_try.exists())
                                {

                                    sb_threshold.append("<a href='");
                                    //sb_threshold.append(f_try.getAbsolutePath());
                                    sb_threshold.append("TFs"+File.separator+f_try.getName());
                                    sb_threshold.append("' target='_blank'><button class=\"button\">"+s.toUpperCase()+"</button>");
                                    sb_threshold.append("</a>");

                                    total_number_tfs.add(s.toUpperCase());
                                    total_number_tfs_hm.add(s.toUpperCase());
                                    total_numbers_tfs_hm_same.add(s.toUpperCase());

                                    HashMap<String,HashSet<String>> current_tf_hm;

                                    if(distinct_tf_hm_diff_same.containsKey(s.toUpperCase()))
                                    {
                                        current_tf_hm=distinct_tf_hm_diff_same.get(s.toUpperCase());
                                    }
                                    else
                                    {
                                        current_tf_hm = new HashMap<>();
                                    }

                                    HashSet<String> current_tf_hm_stage = new HashSet<>();

                                    if(current_tf_hm.containsKey(fileDir_hm.getName()))
                                    {
                                        current_tf_hm_stage = current_tf_hm.get(fileDir_hm.getName());
                                    }
                                    else
                                    {
                                        current_tf_hm_stage= new HashSet<>();
                                    }

                                    current_tf_hm_stage.add("SAME_TPS");
                                    current_tf_hm.put(fileDir_hm.getName(),current_tf_hm_stage);
                                    distinct_tf_hm_diff_same.put(s.toUpperCase(),current_tf_hm);
                                }
                                else
                                {
                                    sb_threshold.append(s.toUpperCase());
                                }

                            }
                            sb_threshold.append("\t\t\t\t</th>\n");
                            count++;
                        }
                        tf_gene_count.put(tf_key.toUpperCase(),tf_key_set);

                        sb_threshold.append("\t\t\t</tr>\n");

                    }
                    br_gc_same.close();

                    sb_threshold.append("\t\t</table>\n");

                    //collapse thing
                    sb_threshold.append("</div>\n");
                    sb_threshold.append("</button>");

                    sb_threshold.append("<h5> A total number of "+total_numbers_tfs_hm_same.size()+" distinct TFs are considered in same time points group. </h5>\n");

                    sb_threshold.append("<h5> A total number of "+total_number_tfs_hm.size()+" distinct TFs are considered in Histone Modification "+fileDir_hm.getName()+ " group. </h5>\n");

                    //collapse thing
                    sb_threshold.append("<button class=\"button_expandable\" id=\"button_"+fileDir_hm.getName()+"_distTFs\" aria-expanded=\"false\" ondblclick=\"expand_collapse('button_"+fileDir_hm.getName()+"_distTFs','table_"+fileDir_hm.getName()+"_distTFs')\">"+fileDir_hm.getName()+": Distinct TFs\n");
                    //tooltip
                    sb_threshold.append("<div class=\"tooltip\"><img src=\".."+File.separator+".."+File.separator+options_intern.folder_out_website_basics+File.separator+options_intern.folder_out_website_basics_website+File.separator+options_intern.folder_out_website_basics_website_images+File.separator+"information.png"+"\" style=\"width:45px;height:40px;\"/>"+
                            "  <span class=\"tooltiptext\">Distinct TFs in a group are the TFs in the union of different and same conditions.</span>\n" +
                            "</div>");
                    sb_threshold.append("<div class='container-buttons' style=\"display: none;background-color: white;color:black;width:100%;\" id=\"table_"+fileDir_hm.getName()+"_distTFs\">\n");
                    sb_threshold.append("\t\t<table style=\"width:100%;font-size:15px;\">\n");
                    for(String tf_distinct : total_number_tfs_hm)
                    {
                        sb_threshold.append("<tr><th><a href='TFs"+File.separator+tf_distinct.toUpperCase()+".html' target='_blank'> <button class=\"button\">"+tf_distinct.toUpperCase()+"</button></a></th></tr>\n");
                    }
                    sb_threshold.append("</table>");

                    //collapse thing
                    sb_threshold.append("</div>\n");
                    sb_threshold.append("</button>\n");



                    sb_threshold.append("    </div>\n" +
                            "\n" +
                            "    <div class=\"w3-third w3-center\">\n" +
                            "      <i class=\"fa fa-anchor w3-padding-64 w3-text-red\"></i>\n");
                    sb_threshold.append("    </div>\n" +
                            "  </div>\n" +
                            "</div>\n");
                }


                sb_threshold.append("<div class=\"container-buttons w3-content\"><h5> A total number of "+total_number_tfs.size()+" distinct TFs are considered in Threshold "+fileDir.getName()+ " group. </h5>\n");

                //collapse thing
                sb_threshold.append("<button class=\"button_expandable\" id=\"button_distTFs\" aria-expanded=\"false\" ondblclick=\"expand_collapse('button_distTFs','table_distTFs')\"> Distinct TFs overall\n");
                //tooltip
                sb_threshold.append("<div class=\"tooltip\"><img src=\".."+File.separator+".."+File.separator+options_intern.folder_out_website_basics+File.separator+options_intern.folder_out_website_basics_website+File.separator+options_intern.folder_out_website_basics_website_images+File.separator+"information.png"+"\" style=\"width:45px;height:40px;\"/>"+
                        "  <span class=\"tooltiptext\">Distinct TFs in a threshold are the TFs in the union of all timepoints and different and same conditions.</span>\n" +
                        "</div>");
                sb_threshold.append("<div class='container-buttons' style=\"display: none;background-color: white;color:black;width:100%;\" id=\"table_distTFs\">\n");
                sb_threshold.append("\t\t<table style=\"width:100%;font-size:15px;\">\n");
                for(String tf_distinct : total_number_tfs)
                {
                    sb_threshold.append("<tr><th><a href='TFs"+File.separator+tf_distinct.toUpperCase()+".html' target='_blank'> <button class=\"button\">"+tf_distinct.toUpperCase()+"</button></a></th></tr>\n");
                }
                sb_threshold.append("</table>");

                //collapse thing
                sb_threshold.append("</div>\n");
                sb_threshold.append("</button>\n");

                sb_threshold.append(write_table_html(Double.parseDouble(fileDir.getName()),"THRESHOLD"));

                sb_threshold.append("</div>\n");
                sb_threshold.append("</div>\n");
                sb_threshold.append("</div>\n");

                sb_threshold.append(html_tail);


                BufferedWriter bw_thresholds = new BufferedWriter(new FileWriter(fileDir.getAbsolutePath()+File.separator+"threshold_"+fileDir.getName()+"_overview.html"));
                bw_thresholds.write(sb_threshold.toString());bw_thresholds.close();

                ArrayList<String> poss_hms_list = new ArrayList(possible_hms);
                ArrayList<String> poss_tfs_list = new ArrayList<>(tfs_to_create_pages);

                File f_output_tf_root = new File(fileDir.getAbsolutePath()+File.separator+options_intern.folder_out_website_htmls_TFs);
                f_output_tf_root.mkdir();

                //create TF pages
                File f_root_target_genes = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_out_target_genes);

                for(String tf : tfs_to_create_pages)
                {
                    File f_output_tf_page = new File(f_output_tf_root.getAbsolutePath()+File.separator+tf.toUpperCase()+".html");
                    if(f_output_tf_page.exists())
                    {
                        //continue;
                    }

                    StringBuilder sb_tf_page = new StringBuilder();
                    sb_tf_page.append(get_header_html("TF",options_intern.analysis_types_regression_coefficient_analysis));
                    sb_tf_page.append(" <script>\n" +
                            " document.title = \"TF: "+tf.toUpperCase()+" TH: "+fileDir.getName()+"\";\n" +
                            " </script>\n");

                    ArrayList<File> files_to_consider = new ArrayList<>();

                    sb_tf_page.append("<div class=\"w3-row-padding w3-padding-64 w3-container\">\n" +
                            "  <div class=\"w3-content\">\n" +
                            "    <div class=\"w3-twothird\">\n" +
                            "      <h1>Transcription Factor: ");
                    sb_tf_page.append(tf.toUpperCase());
                    sb_tf_page.append(" <i>details</i></h1>\n");
                    sb_tf_page.append("<h3>Threshold: "+fileDir.getName().toUpperCase()+"</h3>\n");
                    sb_tf_page.append("<h4><i>Click <a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene="+tf.toUpperCase()+"' target='_blank'><button class=\"button\">"+tf.toUpperCase()+"</button></a> to go to GeneCards</i></h4>\n");



                    sb_tf_page.append("<h4>Gene Count threshold: "+options_intern.plot_cutoff_gcs+"</h4>");
                    sb_tf_page.append("\t\t<table style=\"width:100%\">\n");

                    ArrayList<String> table_header = tf_gene_count.get("TF");
                    sb_tf_page.append("\t\t\t<tr>\n");
                    sb_tf_page.append("\t\t\t\t<th>TF</th>\n");
                    for(String s: table_header)
                    {
                        sb_tf_page.append("\t\t\t\t<th>"+s+"</th>\n");
                    }
                    sb_tf_page.append("\t\t\t</tr>\n");

                    ArrayList<String> table_header_tf = tf_gene_count.get(tf.toUpperCase());
                    sb_tf_page.append("\t\t\t<tr>\n");
                    sb_tf_page.append("\t\t\t\t<th>"+tf.toUpperCase()+"</th>\n");
                    for(String s: table_header_tf)
                    {
                        sb_tf_page.append("\t\t\t\t<th>"+s+"</th>\n");
                    }
                    sb_tf_page.append("\t\t\t</tr>\n");

                    sb_tf_page.append("</table>");

                    for(String hm: possible_hms)
                    {

                        File f_interactive_plots = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_out_website+File.separator+options_intern.folder_out_website_interactive_plots+File.separator+fileDir.getName()+File.separator+hm+File.separator+options_intern.folder_out_website_interactive_plots_tps +File.separator);

                        sb_tf_page.append("<div class=\"w3-row-padding w3-padding-64 w3-container\">\n" +
                                "  <div class=\"w3-content\">\n" +
                                "    <div class=\"w3-twothird\">\n" +
                                "      <h1>");
                        sb_tf_page.append(hm);
                        sb_tf_page.append("</h1>\n");

                        sb_tf_page.append("<h3> Group plots: </h3>\n");
                        sb_tf_page.append("</div>\n");
                        sb_tf_page.append("<button class=\"button_expandable\" style=\"width:1200px\" id=\"button_group_plots_"+hm+"\" aria-expanded=\"false\" ondblclick=\"expand_collapse('button_group_plots_"+hm+"','table_group_plots_"+hm+"')\">"+hm+ " group plots\n");
                        sb_tf_page.append("<div class='container-buttons' style=\"display: none;background-color: white;color:black;width:100%;\" id=\"table_group_plots_"+hm+"\">\n");

                        sb_tf_page.append("\t\t<table style=\"width:100%;font-size:15px;\">\n");

                        for(File fileDir_plot : f_interactive_plots.listFiles())
                        {
                            String[] split_dot = fileDir_plot.getName().split("\\.");
                            String[] split_name = split_dot[0].split("_");

                            sb_tf_page.append("\t  \n" +
                                    "\t <tr><th><h4>"+split_name[1]+ " VS " + split_name[2] +" </h4></th></tr> \n" +
                                    "\t  \n" +
                                    "<tr><th><div class=\"container\" style=\"\n" +
                                    "    width: 980px;\n" +
                                    "    position: relative;\n" +
                                    "    display: block;\n" +
                                    "    overflow-x: scroll;overflow-y: scroll;\">\n" +
                                    "\t  <iframe  id=\"igraph_"+fileDir_plot.getName()+"\" scrolling=\"no\" style=\"border:none;\" seamless=\"seamless\" src=\"");
                            sb_tf_page.append(".."+File.separator+".."+File.separator+".."+File.separator+options_intern.folder_out_website_interactive_plots+File.separator+fileDir.getName()+File.separator+hm+File.separator+options_intern.folder_out_website_interactive_plots_tps +File.separator+fileDir_plot.getName());
                            sb_tf_page.append("\" height=\"400\" width=\"2500\" overflow=\"scroll\"></iframe>\n");
                            sb_tf_page.append("</div></th></tr>\n");
                        }

                        sb_tf_page.append("</table>");


                        sb_tf_page.append("</div>\n</button>\n");

                        sb_tf_page.append("<div class=\"w3-twothird\">\n");




                        //include target genes
                        HashSet<String> already_found_tps = new HashSet<>();

                        File f_consider_different = new File(f_root_target_genes.getAbsolutePath()+File.separator+hm+File.separator+fileDir.getName()+File.separator+"all_data_different");
                        for(File tps_different : f_consider_different.listFiles())
                        {
                            if(tps_different.exists() && !already_found_tps.contains(tps_different.getName()))
                            {
                                files_to_consider.add(tps_different);
                                already_found_tps.add(tps_different.getName());
                            }
                        }

                        File f_consider_same = new File(f_root_target_genes.getAbsolutePath()+File.separator+hm+File.separator+fileDir.getName()+File.separator+"all_data_same");
                        for(File tps_different : f_consider_same.listFiles())
                        {
                            if(tps_different.exists() && !already_found_tps.contains(tps_different.getName()))
                            {
                                files_to_consider.add(tps_different);
                                already_found_tps.add(tps_different.getName());
                            }
                        }
                        Boolean anything_to_write = false;
                        for(File f_considered_no_tf : files_to_consider)
                        {
                            File f_considered = new File(f_considered_no_tf.getAbsolutePath()+File.separator+tf+".csv");
                            if(f_considered.exists())
                            {
                                anything_to_write=true;
                            }
                        }
                        if(anything_to_write)
                        {
                            //tooltip

                            sb_tf_page.append("<h3> Target Genes - top "+options_intern.plot_top_k_genes+" (normalised value):\n");
                            sb_tf_page.append("<div class=\"tooltip\"><img src=\".."+File.separator+".."+File.separator+".."+File.separator+options_intern.folder_out_website_basics+File.separator+options_intern.folder_out_website_basics_website+File.separator+options_intern.folder_out_website_basics_website_images+File.separator+"information.png"+"\" style=\"width:35px;height:30px;\"/>"+
                                    "  <span class=\"tooltiptext\">TargetGenes are harvested from TEPIC output. The higher the normalized affinity value the higher is the ranking.</span>\n" +
                                    "</div></h3>\n");
                            sb_tf_page.append("<p><i>Click on Symbol for GeneCard</i></p>\n");
                        }

                        for(File f_considered_no_tf : files_to_consider)
                        {
                            File f_considered = new File(f_considered_no_tf.getAbsolutePath()+File.separator+tf+".csv");
                            if(f_considered.exists())
                            {
                                sb_tf_page.append("<button class=\"button_expandable\" style=\"width:1200px\" id=\"button_"+hm+"_"+f_considered_no_tf.getName()+"\" aria-expanded=\"false\" ondblclick=\"expand_collapse('button_"+hm+"_"+f_considered_no_tf.getName()+"','table_"+hm+"_"+f_considered_no_tf.getName()+"')\">"+hm+": "+f_considered_no_tf.getName()+"\n");
                                sb_tf_page.append("<div class='container-buttons' style=\"display: none;background-color: white;color:black;width:100%;\" id=\"table_"+hm+"_"+f_considered_no_tf.getName()+"\">\n");


                                sb_tf_page.append("<h4> Point: "+f_considered_no_tf.getName()+" </h4>\n");

                                sb_tf_page.append("\t\t<table style=\"width:100%\">\n");

                                BufferedReader br_point_target_genes = new BufferedReader(new FileReader(f_considered));
                                String line_point_target_getnes ="";
                                while((line_point_target_getnes=br_point_target_genes.readLine())!=null)
                                {
                                    String[] split = line_point_target_getnes.split("\t");

                                    sb_tf_page.append("\t\t\t<tr>\n");

                                    int i = 0;
                                    for(String xx:split)
                                    {
                                        if(xx.equals("NOT_AVAILABLE"))
                                        {
                                            sb_tf_page.append("\t\t\t\t<th>");
                                            sb_tf_page.append("-");
                                            sb_tf_page.append("\t\t\t\t</th>\n");
                                        }
                                        else
                                        {
                                            sb_tf_page.append("\t\t\t\t<th>");
                                            if(i==1&&!xx.equals("SYMBOL"))
                                            {
                                                sb_tf_page.append("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene="+xx.toUpperCase()+"' target='_blank'><button class=\"button\">"+xx.toUpperCase()+"</button></a>");
                                            }
                                            else
                                            {
                                                sb_tf_page.append(xx.toUpperCase());
                                            }
                                            sb_tf_page.append("\t\t\t\t</th>\n");
                                        }
                                        i++;

                                    }
                                    sb_tf_page.append("\t\t\t</tr>\n");
                                }
                                br_point_target_genes.close();

                                sb_tf_page.append("\t\t</table>\n");
                                sb_tf_page.append("</div>\n</button>\n");
                            }
                        }
                    }


                    sb_tf_page.append(html_tail);
                    BufferedWriter bw_tf = new BufferedWriter(new FileWriter(f_output_tf_page));
                    bw_tf.write(sb_tf_page.toString());
                    bw_tf.close();
                }

                File f_output_distribution_analysis = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_out_distribution);
                f_output_distribution_analysis.mkdir();

                File f_output_distribution_analysis_analyzed_tfs = new File(f_output_distribution_analysis.getAbsolutePath()+File.separator+options_intern.folder_out_distribution_analyzed_tfs);
                f_output_distribution_analysis_analyzed_tfs.mkdir();

                BufferedWriter bw_distinct_tf_hm_diff_same = new BufferedWriter(new FileWriter(new File(f_output_distribution_analysis_analyzed_tfs.getAbsolutePath()+File.separator+options_intern.file_suffix_distribution_analysis_analysed_tfs)));

                bw_distinct_tf_hm_diff_same.write("TF\tHM\tSTAGES");
                bw_distinct_tf_hm_diff_same.newLine();

                for(String key_tf : distinct_tf_hm_diff_same.keySet())
                {
                    HashMap<String,HashSet<String>> hm_diff_stages = distinct_tf_hm_diff_same.get(key_tf);

                    for(String key_hm : hm_diff_stages.keySet())
                    {
                        StringBuilder sb = new StringBuilder();
                        sb.append(key_tf);
                        sb.append("\t");
                        sb.append(key_hm);
                        sb.append("\t");

                        int i=0;
                        for(String stages:hm_diff_stages.get(key_hm))
                        {
                            if(i==0)
                            {
                                sb.append(stages);
                            }
                            else
                            {
                                sb.append(";");
                                sb.append(stages);
                            }
                            i++;
                        }

                        bw_distinct_tf_hm_diff_same.write(sb.toString());
                        bw_distinct_tf_hm_diff_same.newLine();

                    }

                }

                bw_distinct_tf_hm_diff_same.close();
            }
        }

        logger.logLine("[WEBSITE] Finished creating overview website.");
    }

    /**
     * get target genes of tfs
     */
    public void get_top_k_target_genes_plots() throws IOException {
        logger.logLine("[PLOTS-TARGET-GENES] Start fetching top "+options_intern.plot_top_k_genes+" target genes for TFs under thresholds.");

        File f_input = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_out_data_plots);
        File f_output = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_out_target_genes);
        f_output.mkdir();

        File f_input_target_genes = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tepic_postprocessing+File.separator+options_intern.folder_name_tepic_postprocessing_output);

        HashMap<String,String> ensg_gene_symbol_map = new HashMap<>();
        BufferedReader br_ensg_gene_symbol = new BufferedReader(new FileReader(new File(options_intern.tepic_ensg_symbol)));
        String line_ensg_gene_symbol = br_ensg_gene_symbol.readLine();
        while((line_ensg_gene_symbol=br_ensg_gene_symbol.readLine())!=null)
        {
            String[] split = line_ensg_gene_symbol.split("\t");
            if(split.length>1)
            {
                ensg_gene_symbol_map.put(split[0],split[1]);
            }
        }
        br_ensg_gene_symbol.close();

        for(File fileDir:f_input.listFiles())
        {
            if(fileDir.isDirectory())
            {
                String hm_name = fileDir.getName();

                File f_output_hm = new File(f_output.getAbsolutePath()+File.separator+hm_name);
                f_output_hm.mkdir();

                for (File fileDir_th : fileDir.listFiles())
                {
                    if(fileDir_th.isDirectory())
                    {
                        String th = fileDir_th.getName();

                        File f_output_hm_th = new File(f_output_hm.getAbsolutePath()+File.separator+th);
                        f_output_hm_th.mkdir();

                        for(File fileDir_th_f : fileDir_th.listFiles())
                        {
                            if(fileDir_th_f.isFile())
                            {
                                String file_name = fileDir_th_f.getName();

                                File f_out_hm_th_file = new File(f_output_hm_th.getAbsolutePath()+File.separator+file_name.split("\\.")[0]);
                                f_out_hm_th_file.mkdir();

                                BufferedReader br_tf = new BufferedReader(new FileReader(fileDir_th_f));
                                String line_tf = br_tf.readLine();
                                String[] split_header = line_tf.split(",");
                                ArrayList<String> timepoints_in_header = new ArrayList<>();
                                for(int i = 1; i < split_header.length;i++)
                                {
                                    String[] split_inner = split_header[i].split("VS");

                                    String[] split_left_side = split_inner[0].split(":");

                                    timepoints_in_header.add(split_left_side[1].trim());
                                    timepoints_in_header.add(split_inner[1].trim());
                                }

                                while((line_tf=br_tf.readLine())!=null)
                                {
                                    String[] split = line_tf.split(",",-1);

                                    ArrayList<String> timepoints_identified = new ArrayList<>();

                                    for(int i =1; i < split.length; i++)
                                    {
                                        if(!split[i].equals(""))
                                        {
                                            timepoints_identified.add(timepoints_in_header.get(i-1));
                                            timepoints_identified.add(timepoints_in_header.get(i));
                                        }
                                    }

                                    HashSet<String> already_done_tps = new HashSet<>();

                                    for(int i = 0; i < timepoints_identified.size(); i+=2)
                                    {
                                        if(timepoints_identified.get(i).equals(timepoints_identified.get(i+1)))
                                        {
                                            continue;
                                        }

                                        if(!already_done_tps.contains(timepoints_identified.get(i)))
                                        {
                                            File f_input_target_genes_hm_group_clash = new File(f_input_target_genes.getAbsolutePath()+File.separator+hm_name+File.separator+timepoints_identified.get(i)+"_"+timepoints_identified.get(i+1));

                                            write_target_genes_of_tf(f_input_target_genes_hm_group_clash,timepoints_identified.get(i),f_out_hm_th_file,split[0],ensg_gene_symbol_map);

                                            already_done_tps.add(timepoints_identified.get(i));
                                        }

                                        if(!already_done_tps.contains(timepoints_identified.get(i+1)))
                                        {
                                            File f_input_target_genes_hm_group_clash = new File(f_input_target_genes.getAbsolutePath()+File.separator+hm_name+File.separator+timepoints_identified.get(i)+"_"+timepoints_identified.get(i+1));

                                            write_target_genes_of_tf(f_input_target_genes_hm_group_clash,timepoints_identified.get(i+1),f_out_hm_th_file,split[0],ensg_gene_symbol_map);

                                            already_done_tps.add(timepoints_identified.get(i+1));
                                        }
                                    }
                                }
                                br_tf.close();
                            }
                        }
                    }
                }
            }
        }

        logger.logLine("[PLOTS-TARGET-GENES] Finished fetching top "+options_intern.plot_top_k_genes+" target genes for TFs under thresholds.");
    }

    /**
     * analyze joined dataframe data for interesting TFs
     */
    public void analyze_plots_data() throws IOException {

        logger.logLine("[PLOT-ANALYSE] Analyse plot data.");

        File input_read_counts = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_deseq2_preprocessing+File.separator+options_intern.folder_name_deseq2_preprocessing_gene_symbols);
        File input_data = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_out_data_plots);

        File output_analysis_root = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_out_analysis_data);
        output_analysis_root.mkdir();

        File output_analysis_tp_level = new File(output_analysis_root.getAbsolutePath()+File.separator+options_intern.folder_out_analysis_data_TP_LEVEL);
        output_analysis_tp_level.mkdir();
        File output_analysis_hm_level = new File(output_analysis_root.getAbsolutePath()+File.separator+options_intern.folder_out_analysis_data_HM_LEVEL);
        output_analysis_hm_level.mkdir();
        File output_analysis_website_overview = new File(output_analysis_root.getAbsolutePath()+File.separator+options_intern.folder_out_analysis_data_WEBSITE_OVERVIEW);
        output_analysis_website_overview.mkdir();

        HashMap<String,HashMap<String,HashMap<String,HashMap<String,Double>>>> cutoff_hm_tf_counts_DIFFERENT = new HashMap<>();
        HashMap<String,HashMap<String,HashMap<String,HashMap<String,Double>>>> cutoff_hm_tf_counts_SAME = new HashMap<>();

        HashMap<String,String> composed_tfs = new HashMap<>();
        HashMap<String,HashSet<String>> composed_tfs_tfs = new HashMap<>();
        BufferedReader br_composed_tfs = new BufferedReader(new FileReader(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tepic_postprocessing+File.separator+options_intern.folder_name_tepic_postprocessing_tfs+File.separator+options_intern.file_suffix_tepic_postprocessing_tfs_tfs));
        String line_composed_tfs = "";
        while((line_composed_tfs=br_composed_tfs.readLine())!=null)
        {
            String[] split = line_composed_tfs.split("\t");

            HashSet temp_set = new HashSet();
            for(int i = 1; i < split.length;i++)
            {
                composed_tfs.put(split[i],split[0]);
                temp_set.add(split[i]);
            }
            composed_tfs_tfs.put(split[0],temp_set);
        }
        br_composed_tfs.close();

        HashMap<String,HashMap<String,Double>> tp_tf_gene_counts = new HashMap<>();

        for(File fileDir:input_read_counts.listFiles())
        {
            String tp_name = fileDir.getName().split("\\.")[0];

            HashMap<String,Double> n_hm = new HashMap<>();
            HashMap<String,Double> composed_tfs_counts = new HashMap<>();

            BufferedReader br = new BufferedReader(new FileReader(fileDir));
            String line = br.readLine();
            while((line=br.readLine())!=null)
            {
                String[] split = line.split("\t");

                if(composed_tfs.containsKey(split[0].toUpperCase()))
                {
                    String key_composed_tf = composed_tfs.get(split[0].toUpperCase());
                    double sum=0;

                    if(composed_tfs_counts.containsKey(key_composed_tf))
                    {
                        sum += composed_tfs_counts.get(key_composed_tf);
                    }

                    sum+=Double.parseDouble(split[2]);
                    composed_tfs_counts.put(key_composed_tf,sum);
                }
                else
                {
                    n_hm.put(split[0].toUpperCase(),Double.parseDouble(split[2]));
                }
            }
            br.close();

            n_hm.putAll(composed_tfs_counts);

            tp_tf_gene_counts.put(tp_name,n_hm);

        }


        HashMap<String,HashMap<String,HashMap<String,Boolean>>> hm_th_tf_found= new HashMap<>();

        for(File fileDir: input_data.listFiles())
        {
            if(fileDir.isDirectory())
            {
                String hm = fileDir.getName();
                File output_hm = new File(output_analysis_tp_level.getAbsolutePath()+File.separator+hm);
                output_hm.mkdir();

                File output_website_overview_hm = new File(output_analysis_website_overview.getAbsolutePath()+File.separator+hm);
                output_website_overview_hm.mkdir();


                for(File fileDirHM_th :fileDir.listFiles())
                {
                    if(fileDirHM_th.isDirectory())
                    {
                        HashMap<String,HashMap<String,Boolean>> th_tf_found;

                        HashMap<String,Boolean> tf_found;

                        if(hm_th_tf_found.containsKey(hm))
                        {
                            th_tf_found=hm_th_tf_found.get(hm);
                        }
                        else
                        {
                            th_tf_found=new HashMap<>();
                        }

                        if(th_tf_found.containsKey(fileDirHM_th.getName()))
                        {
                            tf_found=th_tf_found.get(fileDirHM_th.getName());
                        }
                        else
                        {
                            tf_found = new HashMap<>();

                            for(String s: options_intern.website_interesting_tfs)
                            {
                                tf_found.put(s,false);
                            }
                            th_tf_found.put(fileDirHM_th.getName(),tf_found);
                        }

                        hm_th_tf_found.put(hm,th_tf_found);

                        String th_name = fileDirHM_th.getName();
                        File output_th_tp_level = new File(output_hm.getAbsolutePath()+File.separator+fileDirHM_th.getName());
                        output_th_tp_level.mkdir();

                        File output_website_overview_hm_th = new File(output_website_overview_hm.getAbsolutePath()+File.separator+fileDirHM_th.getName());
                        output_website_overview_hm_th.mkdir();

                        for(File fileDir_HM_th_data : fileDirHM_th.listFiles())
                        {
                            if(fileDir_HM_th_data.isFile())
                            {
                                String name = fileDir_HM_th_data.getName().split("\\.")[0];

                                HashMap<String,HashMap<String,HashMap<String,Double>>> th;
                                HashMap<String,HashMap<String,Double>> hm_intern;

                                if(name.matches(".*different.*"))
                                {
                                    if(cutoff_hm_tf_counts_DIFFERENT.containsKey(fileDirHM_th.getName()))
                                    {
                                        th=cutoff_hm_tf_counts_DIFFERENT.get(fileDirHM_th.getName());
                                        if(th.containsKey(hm))
                                        {
                                            hm_intern=th.get(hm);
                                        }
                                        else
                                        {
                                            hm_intern=new HashMap<>();
                                        }
                                    }
                                    else
                                    {
                                        th=new HashMap<>();
                                        hm_intern=new HashMap<>();
                                    }
                                }
                                else
                                {
                                    if(cutoff_hm_tf_counts_SAME.containsKey(fileDirHM_th.getName()))
                                    {
                                        th=cutoff_hm_tf_counts_DIFFERENT.get(fileDirHM_th.getName());
                                        if(th.containsKey(hm))
                                        {
                                            hm_intern=th.get(hm);
                                        }
                                        else
                                        {
                                            hm_intern=new HashMap<>();
                                        }
                                    }
                                    else
                                    {
                                        th=new HashMap<>();
                                        hm_intern=new HashMap<>();
                                    }
                                }

                                BufferedWriter bw = new BufferedWriter(new FileWriter(output_th_tp_level.getAbsolutePath()+File.separator+fileDir_HM_th_data.getName()));
                                StringBuilder sb = new StringBuilder();
                                sb.append("TF");
                                for(String k: tp_tf_gene_counts.keySet())
                                {
                                    sb.append("\t");
                                    sb.append(k);
                                }
                                bw.write(sb.toString());
                                bw.newLine();

                                BufferedReader br = new BufferedReader(new FileReader(fileDir_HM_th_data));
                                String line = br.readLine();
                                while((line=br.readLine())!=null)
                                {
                                    String[] split = line.split(",",-1);

                                    if(tf_found.containsKey(split[0].toUpperCase()))
                                    {
                                        tf_found.put(split[0].toUpperCase(),true);
                                    }

                                    int count = 0;
                                    for(String s:split)
                                    {
                                        if(!s.equals(""))
                                        {
                                            count++;
                                        }

                                    }

                                    if(count>options_intern.plot_cutoff_tps)
                                    {
                                        HashMap<String,Double> tf_gc = new HashMap<>();
                                        //tf_gc.put(split[0],tp_tf_gene_counts.get(name).get(split[0]));
                                        StringBuilder sb_intern = new StringBuilder();
                                        sb_intern.append(split[0]);
                                        int count_row = 0;

                                        for(String k: tp_tf_gene_counts.keySet())
                                        {
                                            if(tp_tf_gene_counts.get(k).containsKey(split[0].toUpperCase()))
                                            {
                                                count_row+=tp_tf_gene_counts.get(k).get(split[0].toUpperCase());

                                                tf_gc.put(k,tp_tf_gene_counts.get(k).get(split[0].toUpperCase()));
                                                sb_intern.append("\t");
                                                sb_intern.append(tp_tf_gene_counts.get(k).get(split[0].toUpperCase()));
                                            }
                                        }

                                        if(count_row>=options_intern.plot_cutoff_gcs)
                                        {
                                            hm_intern.put(split[0].toUpperCase(),tf_gc);
                                            //WRITE
                                            bw.write(sb_intern.toString());
                                            bw.newLine();
                                        }
                                    }
                                }
                                br.close();
                                bw.close();

                                th.put(hm,hm_intern);

                                //save hashmaps correspondingly
                                if(name.matches(".*different.*"))
                                {
                                    cutoff_hm_tf_counts_DIFFERENT.put(th_name,th);
                                }
                                else
                                {
                                    cutoff_hm_tf_counts_SAME.put(th_name,th);
                                }
                            }
                        }

                        BufferedWriter bw_tf_key_overview = new BufferedWriter(new FileWriter(output_website_overview_hm_th.getAbsolutePath()+File.separator+options_intern.file_suffix_website_analysis_tf_available));
                        bw_tf_key_overview.write("TF\tAVAILABLE");
                        bw_tf_key_overview.newLine();
                        for(String tf_key : tf_found.keySet())
                        {
                            bw_tf_key_overview.write(tf_key+"\t"+tf_found.get(tf_key));
                            bw_tf_key_overview.newLine();
                        }
                        bw_tf_key_overview.close();

                    }
                }
            }
        }

        for(String k : cutoff_hm_tf_counts_DIFFERENT.keySet())
        {
            File f_out = new File(output_analysis_hm_level.getAbsolutePath()+File.separator+k);
            f_out.mkdir();

            HashMap<String,HashMap<String,HashMap<String,Double>>> available_hms = cutoff_hm_tf_counts_DIFFERENT.get(k);

            ArrayList<HashMap<String,HashMap<String,Double>>> tf_lists = new ArrayList<>();

            for(String kk:available_hms.keySet())
            {
                tf_lists.add(available_hms.get(kk));
            }

            HashMap<String,HashMap<String,Double>> results = new HashMap<>();

            for(int i = 0; i < tf_lists.size();i++)
            {
                for(String tf:tf_lists.get(i).keySet())
                {
                    int count_hms = 0;

                    for(int j = 0; j < tf_lists.size(); j++)
                    {
                        HashMap<String,HashMap<String,Double>> list = tf_lists.get(j);
                        if(list.containsKey(tf))
                        {
                            count_hms++;
                        }
                    }

                    if(count_hms>=options_intern.plot_cutoff_hms)
                    {
                        results.put(tf,tf_lists.get(i).get(tf));
                    }
                }
            }


            BufferedWriter bw = new BufferedWriter(new FileWriter(f_out.getAbsolutePath()+File.separator+options_intern.file_suffix_analysis_plot_data_hm_level_different));
            StringBuilder sb = new StringBuilder();
            ArrayList<String> tp_order = new ArrayList<>();
            sb.append("TF");
            for(String tp : tp_tf_gene_counts.keySet())
            {
                sb.append("\t");
                sb.append(tp);
                tp_order.add(tp);
            }
            bw.write(sb.toString());
            bw.newLine();

            for(String tf : results.keySet())
            {
                StringBuilder sb_tf = new StringBuilder();
                sb_tf.append(tf);
                HashMap<String,Double> tf_info = results.get(tf);

                for(String order: tp_order)
                {
                    sb_tf.append("\t");
                    sb_tf.append(tf_info.get(order));
                }
                bw.write(sb_tf.toString());
                bw.newLine();


            }
            bw.close();
        }

        logger.logLine("[PLOT-ANALYSE] Finished analysing plot data.");
    }

    /**
     * creates python Scripts for all timepoints of all histone modifications, it also creates an overview plot of all different group clashes (e.g. P_L but not L1_L10) for a set threshold plots
     */
    public void create_tp_plots() throws Exception {

        logger.logLine("[PLOTS] start creating / running python scripts for plots");

        File folder_input = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_out_put_DYNAMITE);
        File folder_output = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_plots);

        File data_output = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_out_data_plots);
        data_output.mkdir();

        File interactive_plots_output_root = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_out_website);
        interactive_plots_output_root.mkdir();

        File interactive_plots_output_interactive_plots = new File(interactive_plots_output_root.getAbsolutePath()+File.separator+options_intern.folder_out_website_interactive_plots);
        interactive_plots_output_interactive_plots.mkdir();

        folder_output.mkdir();

        for(File fileDirHM:folder_input.listFiles())
        {
            if(fileDirHM.isDirectory())
            {
                File folder_outputHM = new File(folder_output.getAbsolutePath()+File.separator+fileDirHM.getName());
                folder_outputHM.mkdir();

                File folder_out_data_HM = new File(data_output.getAbsolutePath()+File.separator+fileDirHM.getName());
                folder_out_data_HM.mkdir();


                for(double d: options_intern.plot_th_coefficient)
                {
                    File out_th = new File(folder_outputHM.getAbsolutePath()+File.separator+d);
                    out_th.mkdir();

                    File out_data_th = new File(folder_out_data_HM.getAbsolutePath()+File.separator+d);
                    out_data_th.mkdir();

                    File interactive_plots_output_coeff = new File(interactive_plots_output_interactive_plots.getAbsolutePath()+File.separator+d);
                    interactive_plots_output_coeff.mkdir();

                    File interactive_plots_output_coeff_hm = new File(interactive_plots_output_coeff.getAbsolutePath()+File.separator+fileDirHM.getName());
                    interactive_plots_output_coeff_hm.mkdir();

                    File interactive_plots_output_coeff_hm_tps = new File(interactive_plots_output_coeff_hm.getAbsolutePath()+File.separator+options_intern.folder_out_website_interactive_plots_tps);
                    interactive_plots_output_coeff_hm_tps.mkdir();

                    File interactive_plots_output_coeff_hm_overview = new File(interactive_plots_output_coeff_hm.getAbsolutePath()+File.separator+options_intern.folder_out_website_interactive_plots_overview);
                    interactive_plots_output_coeff_hm_overview.mkdir();

                    StringBuilder sb = new StringBuilder();

                    /*WEBSITE_INTERACTIVE_PLOTS*/
                    sb.append("import pip\n" +
                            "\n" +
                            "def import_or_install(package):\n" +
                            "    try:\n" +
                            "        __import__(package)\n" +
                            "    except ImportError:\n" +
                            "        pip.main(['install', package])\n\n");

                    sb.append("import io\n" +
                            "from base64 import b64encode\n" +
                            "import_or_install(\"plotly.express\")\n" +
                            "import plotly.express as px\n" +
                            "import_or_install(\"dash\")\n" +
                            "import_or_install(\"dash_core_components\")\n" +
                            "import_or_install(\"dash_html_components\")\n" +
                            "import dash_core_components as dcc\n" +
                            "import dash_html_components as html\n" +
                            "from dash.dependencies import Input, Output\n" +
                            "import plotly.graph_objs as go\n\n");

                    sb.append("import_or_install(\"pandas\")\n");
                    sb.append("import_or_install(\"seaborn\")\n");
                    sb.append("import_or_install(\"matplotlib.pyplot\")\n");
                    /*WEBSITE_INTERACTIVE_PLOTS*/

                    sb.append("import pandas as pd\n");
                    sb.append("import seaborn as sns\n");
                    sb.append("import matplotlib.pyplot as plt\n");
                    sb.append("sns.set_context(\"notebook\")\n");
                    sb.append("color = \"#A6CEE3\"\n");
                    sb.append("sns.set_context(\"talk\")\n");
                    sb.append("sns.set_style(\"whitegrid\")\n");

                    HashSet<String> th_group_differentpoints = new HashSet<>();
                    HashSet<String> th_group_samepoints = new HashSet<>();

                    for(File fileDirHM_Group : fileDirHM.listFiles())
                    {

                        String[] group_split = fileDirHM_Group.getName().split("_");
                        String groupname1 = group_split[0];
                        String groupname2 = group_split[1];

                        if(groupname1.charAt(0)!=groupname2.charAt(0))
                        {
                            th_group_differentpoints.add(fileDirHM_Group.getName());
                        }
                        else
                        {
                            th_group_samepoints.add(fileDirHM_Group.getName());
                        }
                        sb.append("plt.figure(figsize=(26, 20))\n");
                        sb.append("#create barplot for group: ");
                        sb.append(fileDirHM.getName());
                        sb.append("-");
                        sb.append(fileDirHM_Group.getName());
                        sb.append("\n");
                        sb.append("#Read data\n");
                        sb.append(fileDirHM_Group.getName());
                        sb.append("=pd.read_table('");
                        sb.append(fileDirHM_Group.getAbsolutePath()+File.separator+options_intern.file_suffix_dynamite_output_to_be_plotted);
                        sb.append("').sort_values(['value'], ascending=False)\n");
                        sb.append("# Remove suffix\n");
                        sb.append(fileDirHM_Group.getName());
                        sb.append("['TF'] = ");
                        sb.append(fileDirHM_Group.getName());
                        sb.append("['TF'].str.split('_', expand=True)[0]\n");
                        sb.append("# Sort and filter\n");
                        sb.append(fileDirHM_Group.getName()+"_temp1");
                        sb.append(" = ");
                        sb.append(fileDirHM_Group.getName());
                        sb.append("[");
                        sb.append(fileDirHM_Group.getName());
                        sb.append("['value'] > ");
                        sb.append(d);
                        sb.append("].set_index('TF')\n");
                        sb.append(fileDirHM_Group.getName()+"_temp2");
                        sb.append(" = ");
                        sb.append(fileDirHM_Group.getName());
                        sb.append("[");
                        sb.append(fileDirHM_Group.getName());
                        sb.append("['value'] < -");
                        sb.append(d);
                        sb.append("].set_index('TF')\n");
                        sb.append(fileDirHM_Group.getName());
                        sb.append(" = ");
                        sb.append("pd.concat([");
                        sb.append(fileDirHM_Group.getName()+"_temp1, ");
                        sb.append(fileDirHM_Group.getName()+"_temp2],axis=0)\n");
                        sb.append(fileDirHM_Group.getName());
                        sb.append(".columns = ['");
                        sb.append(fileDirHM.getName());
                        sb.append(": ");
                        sb.append(groupname1);
                        sb.append(" VS ");
                        sb.append(groupname2);
                        sb.append("']\n");
                        sb.append("if('Peak' in ");
                        sb.append(fileDirHM_Group.getName());
                        sb.append(".index):\n");
                        sb.append("    ");
                        sb.append(fileDirHM_Group.getName());
                        sb.append("=pd.DataFrame.drop(");
                        sb.append(fileDirHM_Group.getName());
                        sb.append(",index='Peak')\n");
                        sb.append("if not ");
                        sb.append(fileDirHM_Group.getName());
                        sb.append(".empty:");
                        sb.append("    # Bar Plot\n");
                        sb.append("    time = '");
                        sb.append(fileDirHM.getName());
                        sb.append(": ");
                        sb.append(groupname1);
                        sb.append(" VS ");
                        sb.append(groupname2);
                        sb.append("'\n");
                        sb.append("    ax = sns.barplot(x= ");
                        sb.append(fileDirHM_Group.getName());
                        sb.append(".index, y=time, data=");
                        sb.append(fileDirHM_Group.getName());
                        sb.append(", color=color)\n");
                        sb.append("    ax.set_title(time)\n");
                        sb.append("    ax.set_ylabel('Normalized feature value')\n");
                        sb.append("    ax.set_xlabel('')\n");
                        sb.append("    plt.xticks(rotation=90)\n");
                        sb.append("    plt.tight_layout()\n");
                        sb.append("    plt.savefig(f\"");
                        sb.append(out_th.getAbsolutePath()+File.separator+fileDirHM.getName()+"_"+fileDirHM_Group.getName()+"_threshold_"+d+".png");
                        sb.append("\")\n");
                        sb.append("    plt.clf()\n");
                        /*WEBSITE_INTERACTIVE_PLOTS*/
                        sb.append("    fig_fin = px.bar(");
                        sb.append(fileDirHM_Group.getName());
                        sb.append(",x= ");
                        sb.append(fileDirHM_Group.getName());
                        sb.append(".index, y=time)\n");
                        sb.append("    fig_fin.layout.xaxis.dtick=1\n");
                        sb.append("    fig_fin.write_html(f\"");
                        sb.append(interactive_plots_output_coeff_hm_tps.getAbsolutePath()+File.separator+fileDirHM.getName()+"_"+fileDirHM_Group.getName()+"_threshold_"+d+".html\",include_plotlyjs=\"cdn\")\n");
                        /*WEBSITE_INTERACTIVE_PLOTS*/

                    }
                    sb.append("plt.figure(figsize=(26, 20))\n");
                    sb.append("# Heatmap different stages\n");
                    //Clear Peak values

                    for(String s: th_group_differentpoints)
                    {
                        sb.append("if('Peak' in ");
                        sb.append(s);
                        sb.append(".index):\n");
                        sb.append("    ");
                        sb.append(s);
                        sb.append("=pd.DataFrame.drop(");
                        sb.append(s);
                        sb.append(",index='Peak')\n");
                    }

                    for(String s: th_group_differentpoints)
                    {
                        sb.append(s);
                        sb.append("=");
                        sb.append(s);
                        sb.append("[~");
                        sb.append(s);
                        sb.append(".index.duplicated(keep='first')]\n");
                    }

                    sb.append("join_df_stages = pd.concat([");
                    int c = 0;
                    for(String s : th_group_differentpoints)
                    {
                        if(c==0)
                        {
                            sb.append(s);
                        }
                        else
                        {
                            sb.append(", ");
                            sb.append(s);
                        }
                        c++;
                    }
                    sb.append("], axis=1)\n");

                    sb.append("join_df_stages.to_csv(r'");
                    sb.append(out_data_th.getAbsolutePath()+File.separator);
                    sb.append(options_intern.file_suffix_analysis_plot_data_hm_level_different);
                    sb.append("',index = True, header = True)");
                    sb.append("\n");

                    sb.append("if not ");
                    sb.append("join_df_stages");
                    sb.append(".empty:\n");

                    sb.append("    plot = sns.heatmap(join_df_stages.transpose(), cmap=\"Paired\",  square=True, vmin=1, vmax=1, cbar=False, linewidths=0.5, linecolor='black', xticklabels=True)\n");
                    sb.append("    plt.savefig(\"");
                    sb.append(out_th.getAbsolutePath()+File.separator+fileDirHM.getName()+"_threshold_"+d+"_different_stages.png\")\n");

                    sb.append("    plt.clf()\n");

                    /*WEBSITE_INTERACTIVE_PLOTS*/
                    sb.append("    fig_fin = px.imshow(join_df_stages.transpose())\n");
                    sb.append("    fig_fin.update_layout(plot_bgcolor='rgb(136, 136, 136)')\n");
                    sb.append("    fig_fin.update_xaxes(gridcolor='rgb(136, 136, 136)')\n");
                    sb.append("    fig_fin.layout.xaxis.dtick=1\n");
                    sb.append("    fig_fin.update_yaxes(gridcolor='rgb(136, 136, 136)')\n");
                    sb.append("    fig_fin.layout.height=500\n");
                    sb.append("    fig_fin.layout.width=2500\n");
                    sb.append("    fig_fin.layout.yaxis.dtick=1\n");
                    //sb.append("    fig_fin.update_layout(legend=dict(xanchor=\"right\", x=-1))\n");
                    sb.append("    fig_fin.write_html(f\"");
                    sb.append(interactive_plots_output_coeff_hm_overview.getAbsolutePath()+File.separator+fileDirHM.getName()+"_threshold_"+d+"_different_stages.html\",include_plotlyjs=\"cdn\")\n");
                    /*WEBSITE_INTERACTIVE_PLOTS*/

                    sb.append("plt.figure(figsize=(26, 20))\n");

                    sb.append("# Heatmap same stages\n");
                    for(String s: th_group_samepoints)
                    {
                        sb.append("if('Peak' in ");
                        sb.append(s);
                        sb.append(".index):\n");
                        sb.append("    ");
                        sb.append(s);
                        sb.append("=pd.DataFrame.drop(");
                        sb.append(s);
                        sb.append(",index='Peak')\n");
                    }

                    for(String s: th_group_samepoints)
                    {
                        sb.append(s);
                        sb.append("=");
                        sb.append(s);
                        sb.append("[~");
                        sb.append(s);
                        sb.append(".index.duplicated(keep='first')]\n");
                    }

                    sb.append("join_df_same = pd.concat([");
                    c = 0;
                    for(String s : th_group_samepoints)
                    {
                        if(c==0)
                        {
                            sb.append(s);
                        }
                        else
                        {
                            sb.append(", ");
                            sb.append(s);
                        }
                        c++;
                    }
                    sb.append("], axis=1)\n");

                    sb.append("join_df_same.to_csv(r'");
                    sb.append(out_data_th.getAbsolutePath()+File.separator);
                    sb.append(options_intern.file_suffix_analysis_plot_data_hm_level_same);
                    sb.append("',index = True, header = True)");
                    sb.append("\n");

                    sb.append("if not ");
                    sb.append("join_df_same");
                    sb.append(".empty:\n");

                    sb.append("    plot = sns.heatmap(join_df_same.transpose(), cmap=\"Paired\",  square=True, vmin=1, vmax=1, cbar=False, linewidths=0.5, linecolor='black', xticklabels=True)\n");
                    sb.append("    plt.savefig(\"");
                    sb.append(out_th.getAbsolutePath()+File.separator+fileDirHM.getName()+"_threshold_"+d+"_same_stages.png\")\n");
                    /*WEBSITE_INTERACTIVE_PLOTS*/
                    sb.append("    fig_fin = px.imshow(join_df_same.transpose())\n");
                    sb.append("    fig_fin.update_layout(plot_bgcolor='rgb(136, 136, 136)')\n");
                    sb.append("    fig_fin.update_xaxes(gridcolor='rgb(136, 136, 136)')\n");
                    sb.append("    fig_fin.layout.xaxis.dtick=1\n");
                    sb.append("    fig_fin.update_yaxes(gridcolor='rgb(136, 136, 136)')\n");
                    sb.append("    fig_fin.layout.height=500\n");
                    sb.append("    fig_fin.layout.width=2500\n");
                    //sb.append("    fig_fin.update_layout(legend=dict(yanchor=\"top\", y=0.99, xanchor=\"left\", x=0.01))\n");
                    sb.append("    fig_fin.layout.yaxis.dtick=1\n");
                    sb.append("    fig_fin.write_html(f\"");
                    sb.append(interactive_plots_output_coeff_hm_overview.getAbsolutePath()+File.separator+fileDirHM.getName()+"_threshold_"+d+"_same_stages.html\",include_plotlyjs=\"cdn\")\n");
                    /*WEBSITE_INTERACTIVE_PLOTS*/


                    BufferedWriter bw = new BufferedWriter(new FileWriter(new File(out_th.getAbsolutePath()+File.separator+fileDirHM.getName()+"_"+d+".py")));
                    bw.write(sb.toString());
                    bw.close();

                    String command = "python3 " + out_th.getAbsolutePath()+File.separator+fileDirHM.getName()+"_"+d+".py";

                    logger.logLine("[PLOTS] run plot: " + command);
                    Process child = Runtime.getRuntime().exec(command);
                    int code = child.waitFor();
                    switch (code) {
                        case 0:
                            break;
                        case 1:
                            String message = child.getErrorStream().toString();
                            throw new Exception(message);
                    }
                }
            }
        }
        logger.logLine("[PLOTS] end creating / running python scripts for plots");
    }

    /**
     * run DYNAMITE.R
     */
    public void run_DYNAMITE() throws Exception {
        logger.logLine("[DYNAMITE] start running DYNAMITE with parameters: ");

        String command_base = "Rscript " + options_intern.path_to_COM2POSE+File.separator+options_intern.directory_for_tepic_DYNAMITE+File.separator+"DYNAMITE.R";

        File folder_input;
        folder_input=new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_output_preprocessing_DYNAMITE+File.separator+options_intern.folder_output_preprocessing_DYNAMITE_prepareClass);


        File folder_output = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_out_put_DYNAMITE);
        folder_output.mkdir();

        for(File fileDirHM : folder_input.listFiles())
        {
            if(fileDirHM.isDirectory())
            {
                File folder_outputHM = new File(folder_output.getAbsolutePath()+File.separator+fileDirHM.getName());
                folder_outputHM.mkdir();

                for(File fileDirHM_Group : fileDirHM.listFiles())
                {
                    if(fileDirHM_Group.isDirectory())
                    {
                        File folder_outputHM_Group = new File(folder_outputHM.getAbsolutePath()+File.separator+fileDirHM_Group.getName());
                        folder_outputHM_Group.mkdir();

                        String input_dir = fileDirHM_Group.getAbsolutePath()+File.separator;
                        String output_dir = folder_outputHM_Group.getAbsolutePath()+File.separator;

                        String command_edited = new String(command_base);
                        command_edited += " --dataDir="+input_dir;
                        command_edited += " --outDir="+output_dir;
                        command_edited += " --out_var="+options_intern.dynamite_out_var;
                        command_edited += " --Ofolds="+options_intern.dynamite_Ofolds;
                        command_edited += " --Ifolds="+options_intern.dynamite_Ifolds;
                        command_edited += " --performance="+options_intern.dynamite_performance;
                        command_edited += " --alpha="+options_intern.dynamite_alpha;
                        command_edited += " --cores="+options_intern.dynamite_cores;
                        if(options_intern.dynamite_randomise)
                        {
                            command_edited += " --randomise="+options_intern.dynamite_randomise;
                        }

                        logger.logLine("[DYNAMITE] "+ fileDirHM.getName()+":"+ fileDirHM_Group.getName()+" DYNAMITE.R: " + command_edited);
                        logger.logLine("[DYNAMITE] ... waiting ...");
                        Process child = Runtime.getRuntime().exec(command_edited);
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
            }
        }

        logger.logLine("[DYNAMITE] finished running DYNAMITE");
    }
    /**
     * run integrateData.py and prepareForClassification.R
     */
    public void preprocess_dynamite() throws Exception {

        logger.logLine("[DYNAMITE]: start preprocessing data for DYNAMITE");

        File folder;
        if(options_intern.path_tgen.equals(""))
        {
            folder= new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tepic_postprocessing+File.separator+options_intern.folder_name_tepic_postprocessing_output);
        }
        else
        {
            if(options_intern.tgen_self_regulatory)
            {
                folder= new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tgen+File.separator+options_intern.folder_name_tgen_integrate);
            }
            else
            {
                folder = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tgen+File.separator+options_intern.folder_name_tgen_filter_target_genes);
            }
        }

        File folder_output_pre = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_output_preprocessing_DYNAMITE);
        folder_output_pre.mkdir();

        File folder_output = new File(folder_output_pre.getAbsolutePath()+File.separator+options_intern.folder_output_preprocessing_DYNAMITE_integrateData);
        folder_output.mkdir();

        for(File fileDirHM : folder.listFiles())
        {
            if(fileDirHM.isDirectory())
            {
                File folder_output_hm = new File(folder_output.getAbsolutePath()+File.separator+fileDirHM.getName());
                folder_output_hm.mkdir();

                for(File fileDirGroups : fileDirHM.listFiles())
                {
                    if(fileDirGroups.isDirectory())
                    {
                        File folder_output_group = new File(folder_output_hm.getAbsolutePath()+ File.separator + fileDirGroups.getName());
                        folder_output_group.mkdir();

                        String[] groups_splitted = fileDirGroups.getName().split("_");
                        String groupname1 = groups_splitted[0];
                        String groupname2 = groups_splitted[1];

                        File input_ratios;
                        if(options_intern.path_tgen.equals(""))
                        {
                            input_ratios = new File(fileDirGroups.getAbsolutePath()+File.separator+options_intern.folder_name_tepic_postprocessing_output_ratios+File.separator+options_intern.file_suffix_tepic_postprocessing_output_ratios+fileDirGroups.getName()+".txt");
                        }
                        else
                        {
                            input_ratios = new File(fileDirGroups.getAbsolutePath()+File.separator+options_intern.file_suffix_tepic_postprocessing_output_ratios+groupname1+"_"+groupname2+".txt");
                        }

                        File input_diff_gene_expr = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_deseq2_output+File.separator+fileDirGroups.getName()+options_intern.file_suffix_deseq2_output_DYNAMITE);

                        String command = "python " + options_intern.path_to_COM2POSE+File.separator+options_intern.directory_for_tepic_DYNAMITE+File.separator+"integrateData.py";
                        command += " "+ input_ratios.getAbsolutePath();
                        command += " " + input_diff_gene_expr.getAbsolutePath();
                        command += " " + folder_output_group.getAbsolutePath()+File.separator+options_intern.file_suffix_output_preprocessing_DYNAMITE_integrateData_log2coeff;
                        if(options_intern.dynamite_preprocessing_integrate_data_geneIDs!=0)
                        {
                            command += " " + options_intern.dynamite_preprocessing_integrate_data_geneIDs;
                        }
                        if(options_intern.dynamite_preprocessing_integrate_data_log2fc!=1)
                        {
                            command += " " + options_intern.dynamite_preprocessing_integrate_data_log2fc;
                        }
                        if(!options_intern.dynamite_preprocessing_integrate_data_consider_geneFile.equals(""))
                        {
                            command += " " + options_intern.dynamite_preprocessing_integrate_data_consider_geneFile;
                        }

                        logger.logLine("[DYNAMITE] " +fileDirHM.getName()+":"+ fileDirGroups.getName()+" preprocessing integrateData.py: " + command);
                        Process child = Runtime.getRuntime().exec(command);
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

            }
        }

        File folder_output_classification = new File(folder_output_pre+File.separator+options_intern.folder_output_preprocessing_DYNAMITE_prepareClass);
        folder_output_classification.mkdir();

        for(File fileDirHM:folder_output.listFiles())
        {
            if(fileDirHM.isDirectory())
            {
                File folder_output_classificationHM = new File(folder_output_classification.getAbsolutePath()+File.separator+fileDirHM.getName());
                folder_output_classificationHM.mkdir();

                for(File fileDirGroup:fileDirHM.listFiles())
                {
                    if(fileDirGroup.isDirectory())
                    {
                        File folder_output_classificationHM_Group= new File(folder_output_classificationHM.getAbsolutePath()+File.separator+fileDirGroup.getName());
                        folder_output_classificationHM_Group.mkdir();

                        String command = "Rscript "+ options_intern.path_to_COM2POSE+File.separator+options_intern.directory_for_tepic_DYNAMITE+File.separator+"prepareForClassificiation.R";
                        command += " " + fileDirGroup.getAbsolutePath()+File.separator+options_intern.file_suffix_output_preprocessing_DYNAMITE_integrateData_log2coeff;
                        command += " " + folder_output_classificationHM_Group.getAbsolutePath()+File.separator+options_intern.file_suffix_output_preprocessing_DYNAMITE_prepClass;
                        logger.logLine("[DYNAMITE] " + fileDirGroup.getName()+" preprocessing prepareForClassification.R: " + command);
                        Process child = Runtime.getRuntime().exec(command);
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

            }
        }

        logger.logLine("[DYNAMITE]: finished preprocessing data for DYNAMITE");

    }

    /**
     * integrates self-regulatory TGene TFs into TEPIC data into one affinity file so preproccesing DYNAMITE can use it.
     */
    public void integrate_self_regulatory_tgen() throws IOException {
        logger.logLine("[TGENE-SELF-REGULATORY] Integrate self-regulatoring TFs found in TGene data into TEPIC data");

        File folder_output = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tgen+File.separator+options_intern.folder_name_tgen_integrate);
        folder_output.mkdir();
        //create necessary folder structure -> TEMPLATE: postprocess TEPIC
        File folder_tepic_postprocess = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tgen+File.separator+options_intern.folder_name_tgen_filter_target_genes);

        for(File firDir:folder_tepic_postprocess.listFiles())
        {
            if(firDir.isDirectory())
            {
                File folder_output_hm = new File(folder_output.getAbsolutePath()+File.separator+firDir.getName());
                folder_output_hm.mkdir();
                for(File fileDirHM: firDir.listFiles())
                {
                    if(fileDirHM.isDirectory())
                    {
                        File folder_output_HM_group= new File(folder_output_hm.getAbsolutePath()+File.separator+fileDirHM.getName());
                        folder_output_HM_group.mkdir();
                    }
                }
            }
        }


        //integrate data between TEPIC and TGENE (modifying quotients)
        //based on structure build necessary files
        for(File folderDirHM:folder_output.listFiles())
        {
            if (folderDirHM.isDirectory())
            {
                String hm = folderDirHM.getName();
                for (File folderDirHM_Group : folderDirHM.listFiles())
                {

                    if (folderDirHM_Group.isDirectory())
                    {
                        logger.logLine("[TGENE-SELF-REGULATORY] integrate self-regulatory TFs " + hm +": " + folderDirHM_Group.getName());

                        BufferedWriter bw = new BufferedWriter(new FileWriter(new File(folderDirHM_Group.getAbsolutePath() + File.separator + options_intern.file_suffix_tepic_postprocessing_output_ratios+folderDirHM_Group.getName()+".txt")));

                        File input_data_TGENE = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tgen+File.separator+options_intern.folder_name_tgen_groups+File.separator+hm+File.separator+folderDirHM_Group.getName()+File.separator+options_intern.file_suffic_tgen_output_groups);
                        File input_data_TEPIC = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tepic_postprocessing+File.separator+options_intern.folder_name_tepic_postprocessing_output+File.separator+hm+File.separator+folderDirHM_Group.getName()+File.separator+options_intern.folder_name_tepic_postprocessing_output_ratios+File.separator+options_intern.file_suffix_tepic_postprocessing_output_ratios+folderDirHM_Group.getName()+".txt");

                        //new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_output_preprocessing_DYNAMITE+File.separator+options_intern.folder_output_preprocessing_DYNAMITE_integrateData+File.separator+hm+File.separator+folderDirHM_Group.getName()+File.separator+options_intern.file_suffix_output_preprocessing_DYNAMITE_integrateData_log2coeff);

                        File input_map_position_ensg = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tgen+File.separator+options_intern.folder_name_tgen_preprocessing+File.separator+options_intern.folder_name_tgen_preprocessing_binary_trees+File.separator+options_intern.folder_name_tgen_preprocessing_binary_trees_sorted);

                        HashMap<String,ArrayList<ENSG_ranges_binary_trees>> regions = new HashMap<>();
                        HashMap<String,Boolean> tfs_in_tgene = new HashMap<>();

                        HashMap<String, ENSG_binary_tree> chr_binary_tree = new HashMap<>();

                        for(File chr : input_map_position_ensg.listFiles())
                        {
                            ArrayList<ENSG_ranges_binary_trees> chr_region_list = new ArrayList<>();

                            BufferedReader br = new BufferedReader(new FileReader(chr));
                            String line = br.readLine();
                            while((line=br.readLine())!=null)
                            {
                                String[] split = line.split("\t");

                                ENSG_ranges_binary_trees iu = new ENSG_ranges_binary_trees();
                                iu.chromosome=chr.getName();
                                iu.number = Integer.parseInt(split[0]);
                                iu.left_border=Integer.parseInt(split[1]);
                                iu.right_border=Integer.parseInt(split[2]);
                                iu.ensgs.addAll(Arrays.asList(split[3].split(";")));

                                chr_region_list.add(iu);

                            }
                            br.close();

                            regions.put(chr.getName(),chr_region_list);
                        }

                        for(String chr : regions.keySet())
                        {
                            ArrayList<ENSG_ranges_binary_trees> current_regions = regions.get(chr);

                            ENSG_binary_tree_node root = new ENSG_binary_tree_node(current_regions.get(0),current_regions.get(0).number);
                            ENSG_binary_tree bin_tree = new ENSG_binary_tree(root);

                            for(int i = 1; i < current_regions.size();i++)
                            {
                                bin_tree.add(current_regions.get(i).number,current_regions.get(i));

                            }

                            chr_binary_tree.put(chr.split("\\.")[0],bin_tree);

                        }

                        BufferedReader br_tgene = new BufferedReader(new FileReader(input_data_TGENE));
                        String line_tgene = br_tgene.readLine();

                        HashMap<String,ENSG_ranges_binary_trees> ensgTargetGenes_TFs = new HashMap<>();

                        while((line_tgene= br_tgene.readLine())!=null)
                        {
                            String split[] = line_tgene.split("\t");

                            String chr = split[0];

                            ENSG_ranges_binary_trees iu = new ENSG_ranges_binary_trees();
                            iu.chromosome=chr;
                            iu.left_border= Integer.parseInt(split[1]);
                            iu.right_border=Integer.parseInt(split[2]);
                            iu.ensgs.addAll(Arrays.asList(split[3].split(";")));

                            ENSG_binary_tree tree;

                            if(chr.startsWith("M"))
                            {
                                tree = chr_binary_tree.get("chr"+options_intern.tgen_mt_writing);
                            }
                            else
                            {
                                tree = chr_binary_tree.get("chr"+chr);
                            }

                            if(tree == null)
                            {
                                continue;
                            }

                            ENSG_ranges_binary_trees ensg_match = tree.containsNode(iu);

                            if(ensg_match == null)
                            {
                                continue;
                            }

                            for(String k: ensg_match.ensgs)
                            {
                                if(ensgTargetGenes_TFs.containsKey(k))
                                {
                                    ensgTargetGenes_TFs.get(k).ensgs.addAll(iu.ensgs);
                                }
                                else
                                {
                                    ensgTargetGenes_TFs.put(k,iu);
                                }
                            }

                            for(String k: iu.ensgs)
                            {
                                tfs_in_tgene.put(k.toUpperCase(),false);
                            }
                        }

                        BufferedReader br_tepic = new BufferedReader(new FileReader(input_data_TEPIC));
                        String line_tepic = br_tepic.readLine();

                        bw.write(line_tepic);
                        bw.newLine();

                        HashMap<String,Integer> tf_place = new HashMap<>();

                        ArrayList<String> header = new ArrayList<>(Arrays.asList(line_tepic.split("\t".toUpperCase())));
                        for(int i = 0; i < header.size();i++)
                        {
                            String key = header.get(i).split("_")[0].toUpperCase();
                            header.set(i,key);
                            String[] split_doubles = key.split("::");
                            for(int j = 0; j < split_doubles.length; j++)
                            {
                                tf_place.put(split_doubles[j],i);
                            }
                        }

                        while((line_tepic= br_tepic.readLine())!=null)
                        {
                            String[] split = line_tepic.split("\t");

                            StringBuilder sb = new StringBuilder();
                            sb.append(split[0]);

                            if(ensgTargetGenes_TFs.containsKey(split[0]))
                            {
                                HashSet<String> tgene_pref_tfs =  ensgTargetGenes_TFs.get(split[0]).ensgs;
                                HashSet<Integer> pref_positions = new HashSet<>();

                                for(String k: tgene_pref_tfs)
                                {
                                    pref_positions.add(tf_place.get(k));
                                }

                                for(int i = 1; i < split.length;i++)
                                {
                                    if(pref_positions.contains(i))
                                    {
                                        double score = Double.parseDouble(split[i]);
                                        if(options_intern.tgen_consensus_calc.equals("INCREASE_TGENE_TFS"))
                                        {
                                            score = score + (score * (1-options_intern.tgen_consensus));

                                        }
                                        sb.append("\t");
                                        sb.append(score);
                                    }
                                    else
                                    {
                                        double score = Double.parseDouble(split[i]);
                                        if(options_intern.tgen_consensus_calc.equals("DECREASE_NOT_TGENE_TFs"))
                                        {
                                            score *= (1-options_intern.tgen_consensus);
                                        }
                                        sb.append("\t");
                                        sb.append(score);
                                    }
                                }
                            }
                            else
                            {
                                for(int i = 1; i < split.length;i++)
                                {
                                    double score = Double.parseDouble(split[i]);
                                    if(options_intern.tgen_consensus_calc.equals("DECREASE_NOT_TGENE_TFs"))
                                    {
                                        score *= (1-options_intern.tgen_consensus);
                                    }
                                    sb.append("\t");
                                    sb.append(score);
                                }
                            }

                            bw.write(sb.toString());
                            bw.newLine();
                        }

                        br_tepic.close();
                        bw.close();
                    }
                }
            }
        }

        logger.logLine("[TGENE-SELF-REGULATORY] Finished integrating self-regulatory TFs found in TGene data into TEPIC data");
    }

    /**
     * filter targetgenes of TEPIC using TGENE output
     */
    public void filter_target_genes_tgen() throws IOException {
        logger.logLine("[TGENE] Start filtering target genes.");

        File folder_output = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tgen+File.separator+options_intern.folder_name_tgen_filter_target_genes);
        folder_output.mkdir();
        //create necessary folder structure -> TEMPLATE: postprocess TEPIC
        File folder_tepic_postprocess = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tepic_postprocessing+File.separator+options_intern.folder_name_tepic_postprocessing_output);

        HashMap<String,String> composed_tfs = new HashMap<>();
        BufferedReader br_composed_tfs = new BufferedReader(new FileReader(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tepic_postprocessing+File.separator+options_intern.folder_name_tepic_postprocessing_tfs+File.separator+options_intern.file_suffix_tepic_postprocessing_tfs_tfs));
        String line_composed_tfs = "";
        while((line_composed_tfs=br_composed_tfs.readLine())!=null)
        {
            String[] split = line_composed_tfs.split("\t");

            for(int i = 1; i < split.length;i++)
            {
                composed_tfs.put(split[i],split[0]);
            }
        }
        br_composed_tfs.close();

        for(File firDir:folder_tepic_postprocess.listFiles())
        {
            if(firDir.isDirectory())
            {
                File folder_output_hm = new File(folder_output.getAbsolutePath()+File.separator+firDir.getName());
                folder_output_hm.mkdir();
                for(File fileDirHM: firDir.listFiles())
                {
                    if(fileDirHM.isDirectory())
                    {
                        File folder_output_HM_group= new File(folder_output_hm.getAbsolutePath()+File.separator+fileDirHM.getName());
                        folder_output_HM_group.mkdir();
                    }
                }
            }
        }

        HashMap<String,String> ensg_gene_symbol_map = new HashMap<>();
        BufferedReader br_ensg_gene_symbol = new BufferedReader(new FileReader(new File(options_intern.tepic_ensg_symbol)));
        String line_ensg_gene_symbol = br_ensg_gene_symbol.readLine();
        while((line_ensg_gene_symbol=br_ensg_gene_symbol.readLine())!=null)
        {
            String[] split = line_ensg_gene_symbol.split("\t");
            if(split.length>1)
            {
                if(composed_tfs.containsKey(split[0]))
                {
                    ensg_gene_symbol_map.put(split[1].toUpperCase(),composed_tfs.get(split[0]));
                }
                else
                {
                    ensg_gene_symbol_map.put(split[1].toUpperCase(),split[0]);
                }
            }
        }
        br_ensg_gene_symbol.close();

        for(File folderDirHM:folder_output.listFiles())
        {
            if (folderDirHM.isDirectory())
            {
                String hm = folderDirHM.getName();
                for (File folderDirHM_Group : folderDirHM.listFiles())
                {
                    if (folderDirHM_Group.isDirectory())
                    {
                        String group_name = folderDirHM_Group.getName();

                        logger.logLine("[TGENE] filter target genes " + hm +": " + group_name);

                        File input_data_TGENE = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tgen+File.separator+options_intern.folder_name_tgen_groups+File.separator+hm+File.separator+folderDirHM_Group.getName()+File.separator+options_intern.file_suffic_tgen_output_groups);
                        File input_data_TEPIC = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tepic_postprocessing+File.separator+options_intern.folder_name_tepic_postprocessing_output+File.separator+hm+File.separator+folderDirHM_Group.getName()+File.separator+options_intern.folder_name_tepic_postprocessing_output_ratios+File.separator+options_intern.file_suffix_tepic_postprocessing_output_ratios+folderDirHM_Group.getName()+".txt");

                        HashSet<String> available_target_genes = new HashSet<>();

                        BufferedReader br_tgene = new BufferedReader(new FileReader(input_data_TGENE));
                        String line_tgene = br_tgene.readLine();
                        while((line_tgene = br_tgene.readLine())!=null)
                        {
                            String[] split = line_tgene.split("\t");

                            String[] split_genes = split[3].split(";");

                            for(String s: split_genes)
                            {
                                if(ensg_gene_symbol_map.containsKey(s.toUpperCase()))
                                {
                                    available_target_genes.add(ensg_gene_symbol_map.get(s.toUpperCase()));
                                }
                            }
                        }
                        br_tgene.close();

                        BufferedWriter bw_tepic = new BufferedWriter(new FileWriter(folderDirHM_Group.getAbsolutePath()+File.separator+input_data_TEPIC.getName()));
                        BufferedReader br_tepic = new BufferedReader(new FileReader(input_data_TEPIC));
                        String line_tepic = br_tepic.readLine();
                        bw_tepic.write(line_tepic);
                        bw_tepic.newLine();

                        while((line_tepic=br_tepic.readLine())!=null)
                        {
                            String[] split = line_tepic.split("\t");
                            if(available_target_genes.contains(split[0]))
                            {
                                bw_tepic.write(line_tepic);
                                bw_tepic.newLine();
                            }
                        }

                        br_tepic.close();
                        bw_tepic.close();
                    }
                }
            }
        }
        logger.logLine("[TGENE] Finished filtering target genes.");
    }

    /**
     * creates TGENE group clash so it can be merged into TEPIC
     */
    public void create_tgen_groups() throws IOException {
        logger.logLine("[TGENE] Create TGene groupe data");

        File folder_input = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tgen+File.separator+options_intern.folder_name_tgen_merged);
        File folder_output = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tgen+File.separator+options_intern.folder_name_tgen_groups);
        folder_output.mkdir();

        //create necessary folder structure -> TEMPLATE: postprocess TEPIC
        File folder_tepic_postprocess = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tepic_postprocessing+File.separator+options_intern.folder_name_tepic_postprocessing_output);
        for(File firDir:folder_tepic_postprocess.listFiles())
        {
            if(firDir.isDirectory())
            {
                File folder_output_hm = new File(folder_output.getAbsolutePath()+File.separator+firDir.getName());
                folder_output_hm.mkdir();
                for(File fileDirHM: firDir.listFiles())
                {
                    if(fileDirHM.isDirectory())
                    {
                        File folder_output_HM_group= new File(folder_output_hm.getAbsolutePath()+File.separator+fileDirHM.getName());
                        folder_output_HM_group.mkdir();
                    }
                }
            }
        }

        //based on structure build necessary files
        for(File folderDirHM:folder_output.listFiles())
        {
            if(folderDirHM.isDirectory())
            {
                String hm = folderDirHM.getName();
                for(File folderDirHM_Group : folderDirHM.listFiles())
                {

                    if(folderDirHM_Group.isDirectory())
                    {
                        BufferedWriter bw = new BufferedWriter(new FileWriter(new File(folderDirHM_Group.getAbsolutePath()+File.separator+options_intern.file_suffic_tgen_output_groups)));
                        bw.write("CHR\tLEFT_BORDER\tRIGHT_BORDER\tENSGS");
                        bw.newLine();

                        String[] split_name = folderDirHM_Group.getName().split("_");
                        String group1 = split_name[0];
                        String group2 = split_name[1];

                        File input_group1 = new File(folder_input.getAbsolutePath()+File.separator+group1+File.separator+hm+File.separator+hm+"_"+group1+".txt");
                        File input_group2 = new File(folder_input.getAbsolutePath()+File.separator+group2+File.separator+hm+File.separator+hm+"_"+group2+".txt");

                        HashMap<String,ArrayList<ENSG_ranges_binary_trees>> chr_unmerged_ius = new HashMap<>();

                        BufferedReader br_group1 = new BufferedReader(new FileReader(input_group1));
                        String line_group1 = br_group1.readLine();
                        while((line_group1=br_group1.readLine())!=null)
                        {
                            String[] split_group1 = line_group1.split("\t");
                            String chr = split_group1[0];
                            int left_border = Integer.parseInt(split_group1[1]);
                            int right_border = Integer.parseInt(split_group1[2]);
                            String ensg = split_group1[3].toUpperCase();

                            ArrayList<ENSG_ranges_binary_trees> current_arr_chr;
                            if(chr_unmerged_ius.containsKey(chr))
                            {
                                current_arr_chr = chr_unmerged_ius.get(chr);
                            }
                            else
                            {
                                current_arr_chr = new ArrayList<>();
                            }

                            ENSG_ranges_binary_trees iu = new ENSG_ranges_binary_trees();
                            iu.chromosome=chr;
                            iu.ensgs.add(ensg);
                            iu.left_border=left_border;
                            iu.right_border=right_border;

                            current_arr_chr.add(iu);
                            chr_unmerged_ius.put(chr,current_arr_chr);
                        }
                        br_group1.close();

                        BufferedReader br_group2 = new BufferedReader(new FileReader(input_group2));
                        String line_group2 = br_group2.readLine();
                        while((line_group2=br_group2.readLine())!=null)
                        {
                            String[] split_group2 = line_group2.split("\t");
                            String chr = split_group2[0];
                            int left_border = Integer.parseInt(split_group2[1]);
                            int right_border = Integer.parseInt(split_group2[2]);
                            String ensg = split_group2[3].toUpperCase();

                            ArrayList<ENSG_ranges_binary_trees> current_arr_chr;
                            if(chr_unmerged_ius.containsKey(chr))
                            {
                                current_arr_chr = chr_unmerged_ius.get(chr);
                            }
                            else
                            {
                                current_arr_chr = new ArrayList<>();
                            }

                            ENSG_ranges_binary_trees iu = new ENSG_ranges_binary_trees();
                            iu.chromosome=chr;
                            iu.ensgs.add(ensg);
                            iu.left_border=left_border;
                            iu.right_border=right_border;

                            current_arr_chr.add(iu);
                            chr_unmerged_ius.put(chr,current_arr_chr);
                        }
                        br_group2.close();

                        //sort all arrays
                        for(String s : chr_unmerged_ius.keySet())
                        {
                            ArrayList<ENSG_ranges_binary_trees> x = chr_unmerged_ius.get(s);
                            Collections.sort(x);
                            chr_unmerged_ius.put(s,x);
                        }


                        //remove duplicates or integrate ensgs if in overlap
                        for(String chr:chr_unmerged_ius.keySet())
                        {
                            ArrayList<ENSG_ranges_binary_trees> unmerged = chr_unmerged_ius.get(chr);

                            ArrayList<ENSG_ranges_binary_trees> unmerged_temp = new ArrayList<>();

                            for(int i = 1; i < unmerged.size();i++)
                            {
                                ENSG_ranges_binary_trees iu_before = unmerged.get(i-1);
                                ENSG_ranges_binary_trees iu_after = unmerged.get(i);

                                if(iu_before.isTheSame(iu_after))
                                {
                                    unmerged_temp.add(iu_before);
                                    i++;
                                }
                                else
                                {
                                    unmerged_temp.add(iu_before);
                                }
                            }

                            unmerged = new ArrayList<>(unmerged_temp);
                            unmerged_temp.clear();

                            //check for same_ranges but not same ENSGs and merge them
                            for(int i = 1; i < unmerged.size();i++)
                            {
                                ENSG_ranges_binary_trees iu_before = unmerged.get(i-1);
                                ENSG_ranges_binary_trees iu_after = unmerged.get(i);

                                if(iu_before.isSameRange(iu_after))
                                {
                                    for(String s : iu_after.ensgs)
                                    {
                                        iu_before.ensgs.add(s);
                                    }
                                    unmerged_temp.add(iu_before);

                                    i++;
                                }
                                else
                                {
                                    unmerged_temp.add(iu_before);
                                }
                            }

                            unmerged=new ArrayList<>(unmerged_temp);
                            unmerged_temp.clear();

                            /*
                            boolean something_changed = true;
                            while(something_changed)
                            {
                                something_changed=false;
                                //check for overlaps and if so join them
                                for(int i = 1; i < unmerged.size();i++)
                                {
                                    ENSG_ranges_binary_trees iu_before = unmerged.get(i - 1);
                                    ENSG_ranges_binary_trees iu_after = unmerged.get(i);

                                    if(iu_before.isOverLap(iu_after))
                                    {
                                        iu_before.right_border=iu_after.right_border;
                                        for(String s : iu_after.ensgs)
                                        {
                                            iu_before.ensgs.add(s);
                                        }
                                        something_changed=true;
                                        i++;
                                    }
                                    else
                                    {
                                        unmerged_temp.add(iu_before);
                                    }
                                }
                                unmerged=new ArrayList<>(unmerged_temp);
                                unmerged_temp.clear();
                            }*/


                            for(int i = 0; i < unmerged.size(); i++)
                            {
                                bw.write(unmerged.get(i).toString());
                                bw.newLine();
                            }
                        }
                        bw.close();
                    }
                }
            }
        }
        logger.logLine("[TGENE] Finished creating TGene group data");
    }

    /**
     * merge different samples of one timepoint of TGENE
     */
    public void merge_tgen() throws IOException {
        logger.logLine("[TGENE] Postprocessing: Merge TGene results");

        File folder_input = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tgen+File.separator+options_intern.folder_name_tgen_output);
        File folder_output = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tgen+File.separator+options_intern.folder_name_tgen_merged);
        folder_output.mkdir();

        for(File fileDirTP : folder_input.listFiles())
        {
            if(fileDirTP.isDirectory())
            {
                File output_tp = new File(folder_output.getAbsolutePath()+File.separator+fileDirTP.getName());
                output_tp.mkdir();

                for(File fileDirTP_HM : fileDirTP.listFiles())
                {
                    if(fileDirTP_HM.isDirectory())
                    {
                        File output_tp_hm = new File(output_tp.getAbsolutePath()+File.separator+fileDirTP_HM.getName());
                        output_tp_hm.mkdir();

                        String name_output = fileDirTP_HM.getName() + "_"+ fileDirTP.getName();


                        HashMap<String,ArrayList<ENSG_ranges_binary_trees>> tfs_to_regions = new HashMap<>();

                        for(File fileDirTP_HM_Samples:fileDirTP_HM.listFiles())
                        {
                            if(fileDirTP_HM_Samples.isDirectory())
                            {
                                File input_data = new File(fileDirTP_HM_Samples.getAbsolutePath()+File.separator+options_intern.file_suffix_tgen_output);

                                BufferedReader br = new BufferedReader(new FileReader(input_data));
                                String line = br.readLine();
                                while((line=br.readLine())!=null)
                                {
                                    if(line.startsWith("#") || line.equals(""))
                                    {
                                        continue;
                                    }

                                    String[] split = line.split("\t");

                                    ENSG_ranges_binary_trees iu = new ENSG_ranges_binary_trees();
                                    iu.ensgs.add(split[1]);

                                    String[] split_chr = split[6].split(":");
                                    String[] split_borders = split_chr[1].split("-");

                                    iu.chromosome = split_chr[0];
                                    iu.left_border=Integer.parseInt(split_borders[0]);
                                    iu.right_border=Integer.parseInt(split_borders[1]);

                                    if(tfs_to_regions.containsKey(split_chr[0]))
                                    {
                                        ArrayList<ENSG_ranges_binary_trees> x = tfs_to_regions.get(split_chr[0]);
                                        x.add(iu);
                                        tfs_to_regions.put(split_chr[0],x);

                                    }
                                    else
                                    {
                                        ArrayList<ENSG_ranges_binary_trees> x = new ArrayList<>();
                                        x.add(iu);
                                        tfs_to_regions.put(split_chr[0],x);
                                    }
                                }
                            }
                        }

                        for(String chr: tfs_to_regions.keySet())
                        {
                            Collections.sort(tfs_to_regions.get(chr));
                        }

                        for(String chr: tfs_to_regions.keySet())
                        {
                            ArrayList<ENSG_ranges_binary_trees> tf_reg = tfs_to_regions.get(chr);

                            ArrayList<ENSG_ranges_binary_trees> tf_temp = new ArrayList<>();
                            for(int i = 1; i < tf_reg.size(); i++)
                            {
                                //delete all duplicates
                                ENSG_ranges_binary_trees iu_before = tf_reg.get(i-1);
                                ENSG_ranges_binary_trees iu_now = tf_reg.get(i);

                                if(iu_before.isTheSame(iu_now))
                                {
                                    int temp_i = i+1;
                                    boolean same = true;
                                    int count_same = 0;
                                    while(temp_i < tf_reg.size()&&same)
                                    {
                                        same=false;

                                        ENSG_ranges_binary_trees iu_follow = tf_reg.get(temp_i);

                                        if(iu_follow.isTheSame(iu_before))
                                        {
                                            same=true;
                                            count_same++;
                                        }

                                        temp_i++;
                                    }
                                    i=temp_i+1;


                                    tf_temp.add(iu_before);
                                }
                                else
                                {
                                    tf_temp.add(iu_before);
                                }
                            }
                            tfs_to_regions.put(chr,tf_temp);
                        }

                        //CHECK TPM if TPM filter is set!
                        if(options_intern.tepic_tpm_cutoff>0)
                        {
                                /*command_tail += " -T " + options_intern.tepic_tpm_cutoff;
                                command_tail += " -E " + options_intern.tepic_ensg_symbol;
                                command_tail += " -A " + options_intern.deseq2_input_gene_id;
                                String n_dir = options_intern.com2pose_working_directory+File.separator+ options_intern.folder_name_deseq2_preprocessing+File.separator+options_intern.folder_name_deseq2_preprocessing_single+File.separator+dirGroup.getName()+options_intern.file_suffix_deseq2_preprocessing_meanCounts;
                                command_tail_sample += " -G " + n_dir;

                                [-G input genes count file, if set, default TPM is 1]\n
                                [-T set (T)ranscripts (P)er (M)illion cutoff must be in float form (e.g. 1.0)]
                                [-A input gene annotation desq2 file, required for TPM filter]
                                [-E input ensg to gene symbol file, required for TPM filter*/

                            logger.logLine("[TGENE] TPM filter is set. Filtering TGENE results: " + fileDirTP_HM.getName() + " - " + fileDirTP.getName()+".");

                            double tpm_cutoff = options_intern.tepic_tpm_cutoff;
                            File f_ensg_map_symbol = new File(options_intern.tepic_ensg_symbol);
                            File f_gene_annot_deseq2 = new File(options_intern.deseq2_input_gene_id);
                            File f_gene_count_group1 = new File(options_intern.com2pose_working_directory+File.separator+ options_intern.folder_name_deseq2_preprocessing+File.separator+options_intern.folder_name_deseq2_preprocessing_single+File.separator+fileDirTP.getName()+options_intern.file_suffix_deseq2_preprocessing_meanCounts);
                            File f_ref_genome = new File(options_intern.tepic_gene_annot);


                            ArrayList<Integer> gene_counts = new ArrayList<>();
                            BufferedReader br_gene_counts_1 = new BufferedReader(new FileReader(f_gene_count_group1));
                            String line_gene_count_1 = br_gene_counts_1.readLine();
                            while ((line_gene_count_1= br_gene_counts_1.readLine())!=null)
                            {
                                gene_counts.add(Integer.parseInt(line_gene_count_1));
                            }
                            br_gene_counts_1.close();


                            ArrayList<String> ensg_numbers = new ArrayList<>();
                            BufferedReader br_ensg_numbers = new BufferedReader(new FileReader(f_gene_annot_deseq2));
                            String line_ensg_numbers = br_ensg_numbers.readLine();
                            while((line_ensg_numbers= br_ensg_numbers.readLine())!=null)
                            {
                                ensg_numbers.add(line_ensg_numbers.toUpperCase());
                            }
                            br_ensg_numbers.close();

                            int total_number_samples=0;
                            HashMap<String,Integer> ensg_counts = new HashMap<>();
                            for(int i = 0; i < gene_counts.size();i++)
                            {
                                ensg_counts.put(ensg_numbers.get(i).toUpperCase(),gene_counts.get(i));
                                total_number_samples++;
                            }

                            HashMap<String,String> ensg_symb = new HashMap<>();

                            BufferedReader br_ensg_symb = new BufferedReader(new FileReader(f_ensg_map_symbol));
                            String line_ensg_symb = br_ensg_symb.readLine();
                            while((line_ensg_symb=br_ensg_symb.readLine())!=null)
                            {
                                String[] split = line_ensg_symb.split("\t");
                                if(split.length>1)
                                {
                                    ensg_symb.put(split[1].toUpperCase(),split[0]);
                                }
                            }
                            br_ensg_symb.close();

                            HashMap<String,Integer> gene_symbol_lengths = new HashMap<>();
                            BufferedReader br_gene_symbol_lengths = new BufferedReader(new FileReader(f_ref_genome));
                            String line_gene_symbol_lengths = "";
                            while((line_gene_symbol_lengths=br_gene_symbol_lengths.readLine())!=null)
                            {
                                if(line_gene_symbol_lengths.startsWith("#"))
                                {
                                    continue;
                                }
                                String[] split = line_gene_symbol_lengths.split("\t");
                                if(split[2].equals("transcript"))
                                {
                                    int start=Integer.parseInt(split[3]);
                                    int end=Integer.parseInt(split[4]);
                                    int diff=end-start+1;
                                    String[] split_further = split[8].split(" ");
                                    String ensg_name=split_further[1].replace("\"","");
                                    ensg_name=ensg_name.replace(";","");
                                    String[] ensg_name_split=ensg_name.split("\\.");
                                    gene_symbol_lengths.put(ensg_name_split[0],diff);
                                }
                            }
                            br_gene_symbol_lengths.close();

                            double all_rpk = 0.0;

                            for(String ec : ensg_counts.keySet())
                            {
                                if(!ensg_counts.containsKey(ec))
                                {
                                    continue;
                                }

                                int ec_count = ensg_counts.get(ec);
                                int ec_lengths = 0;
                                if(!gene_symbol_lengths.containsKey(ec))
                                {
                                    continue;
                                }
                                ec_lengths=gene_symbol_lengths.get(ec);
                                if(ec_lengths==0)
                                {
                                    continue;
                                }
                                double norm_ec_lengths= ec_lengths/(1000*1.0);
                                double current_rpk = ec_count/(norm_ec_lengths*1.0);
                                all_rpk+=current_rpk;
                            }

                            double 	scaling_factor=all_rpk/(1000000*1.0);

                            for(String chr: tfs_to_regions.keySet())
                            {
                                ArrayList<ENSG_ranges_binary_trees> tf_reg = tfs_to_regions.get(chr);
                                ArrayList<ENSG_ranges_binary_trees> tf_reg_temp = new ArrayList<>();

                                for(ENSG_ranges_binary_trees iu : tf_reg)
                                {
                                    boolean should_write = false;

                                    //check if we want that TF
                                    for(String ec: iu.ensgs)
                                    {
                                        ec = ec.toUpperCase();

                                        if(ensg_symb.containsKey(ec))
                                        {
                                            ec = ensg_symb.get(ec);
                                        }
                                        else
                                        {
                                            continue;
                                        }

                                        if(!ensg_counts.containsKey(ec))
                                        {
                                            continue;
                                        }

                                        int ec_count = ensg_counts.get(ec);
                                        int ec_lengths = 0;
                                        if(!gene_symbol_lengths.containsKey(ec))
                                        {
                                            continue;
                                        }
                                        ec_lengths=gene_symbol_lengths.get(ec);
                                        if(ec_lengths==0)
                                        {
                                            continue;
                                        }
                                        double norm_ec_lengths= ec_lengths/(1000*1.0);
                                        double current_rpk = ec_count/(norm_ec_lengths*1.0);

                                        if(current_rpk>=scaling_factor)
                                        {
                                            should_write=true;
                                        }
                                    }

                                    if(should_write)
                                    {
                                        tf_reg_temp.add(iu);
                                    }
                                }

                                tfs_to_regions.put(chr,tf_reg_temp);
                            }
                        }


                        BufferedWriter bw = new BufferedWriter(new FileWriter(new File(output_tp_hm.getAbsolutePath()+File.separator+name_output+".txt")));
                        bw.write("CHR\tLEFT_BORDER\tRIGHT_BORDER\tENSG");
                        bw.newLine();
                        for(String chr: tfs_to_regions.keySet())
                        {
                            ArrayList<ENSG_ranges_binary_trees> tf_reg = tfs_to_regions.get(chr);
                            for(ENSG_ranges_binary_trees iu : tf_reg)
                            {
                                bw.write(iu.toString());
                                bw.newLine();
                            }

                        }
                        bw.close();
                    }
                }
            }
        }
        logger.logLine("[TGENE] Finished merging TGene results");
    }

    /**
     * run tgene for all samples
     */
    public void run_tgen() throws Exception {
        logger.logLine("[TGENE] Run TGene");

        String command_base =options_intern.path_tgen+File.separator +"bin"+File.separator+"tgene";

        File output_tgen = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tgen+File.separator+options_intern.folder_name_tgen_output);
        output_tgen.mkdir();

        File folder = new File(options_intern.tepic_input_directory);
        for(File fileDirTP: folder.listFiles())
        {
            if(fileDirTP.isDirectory())
            {
                File output_tgen_tp = new File(output_tgen.getAbsolutePath()+File.separator+fileDirTP.getName());
                output_tgen_tp.mkdir();

                for(File fileDirTP_HM : fileDirTP.listFiles())
                {
                    if(fileDirTP_HM.isDirectory())
                    {
                        File output_tgen_tp_hm = new File(output_tgen_tp.getAbsolutePath()+File.separator+fileDirTP_HM.getName());
                        output_tgen_tp_hm.mkdir();

                        for(File fileDir_data : fileDirTP_HM.listFiles())
                        {
                            if(!fileDir_data.isDirectory())
                            {
                                String[] name_split = fileDir_data.getName().split("\\.");

                                File output_tgen_tp_hm_sample = new File(output_tgen_tp_hm+File.separator+name_split[0]);
                                output_tgen_tp_hm_sample.mkdir();

                                String command_execute = new String(command_base);
                                command_execute += " " + fileDir_data.getAbsolutePath();

                                File gtf_dir = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tgen+File.separator+options_intern.folder_name_tgen_preprocessing+File.separator+options_intern.folder_name_tgen_preprocessing_gtf);
                                File gtf = new File("");

                                for(File gtf_dirDir : gtf_dir.listFiles())
                                {
                                    if(gtf_dirDir.getName().matches(".*"+options_intern.file_suffix_tgen_preprocess_gtf))
                                    {
                                        gtf=gtf_dirDir;
                                        break;
                                    }
                                }

                                if(gtf.getName().equals(""))
                                {
                                    logger.logLine("[TGENE] could not find preprocessed GTF.");
                                    System.exit(1);
                                }

                                command_execute += " " + gtf.getAbsolutePath();

                                command_execute += " -oc " + output_tgen_tp_hm_sample;

                                if(options_intern.tgen_no_closest_locus)
                                {
                                    command_execute += " --no-closest-locus";
                                }
                                if(options_intern.tgen_no_closest_tss)
                                {
                                    command_execute += " --no-closest-tss";
                                }

                                command_execute+= " --max-link-distances " + options_intern.tgen_max_link_distances;
                                command_execute+= " --max-pvalue " + options_intern.tgen_pvalue;

                                //now execute TGENE:
                                logger.logLine("[TGENE] execute TGENE with command line: " + command_execute);
                                Process child = Runtime.getRuntime().exec(command_execute);
                                int code = child.waitFor();
                                switch (code){
                                    case 0:
                                        break;
                                    case 1:
                                        String message = child.getErrorStream().toString();
                                        logger.logLine("[TGENE] please check chromosome names in peak file and gtf file - they need to be named the same!");
                                        throw new Exception(message);
                                }
                            }
                        }
                    }
                }
            }
        }
        logger.logLine("[TGENE] Finished TGene");
    }

    /**
     * preprocess data for tgen run
     */
    public void preprocess_tgen() throws IOException {
        logger.logLine("[TGENE] Consensus approach with TGen is used. Preprocessing ...");

        //create necessary folders for preprocessing
        File f_TGEN = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tgen);
        f_TGEN.mkdir();
        File f_TGEN_preprocess = new File(f_TGEN.getAbsolutePath()+File.separator+options_intern.folder_name_tgen_preprocessing);
        f_TGEN_preprocess.mkdir();
        File f_TGEN_preprocess_gtf = new File(f_TGEN_preprocess.getAbsolutePath()+File.separator+options_intern.folder_name_tgen_preprocessing_gtf);
        f_TGEN_preprocess_gtf.mkdir();
        File f_TGEN_preprocess_binary_trees = new File(f_TGEN_preprocess.getAbsolutePath()+File.separator+options_intern.folder_name_tgen_preprocessing_binary_trees);
        f_TGEN_preprocess_binary_trees.mkdir();
        File f_TGEN_preprocess_binary_trees_unmerged = new File(f_TGEN_preprocess_binary_trees.getAbsolutePath()+File.separator+options_intern.folder_name_tgen_preprocessing_binary_trees_unmerged);
        f_TGEN_preprocess_binary_trees_unmerged.mkdir();
        File f_TGEN_preprocess_binary_trees_merged = new File(f_TGEN_preprocess_binary_trees.getAbsolutePath()+File.separator+options_intern.folder_name_tgen_preprocessing_binary_trees_merged);
        f_TGEN_preprocess_binary_trees_merged.mkdir();
        File f_TGEN_preprocess_binary_trees_sorted = new File(f_TGEN_preprocess_binary_trees.getAbsolutePath()+File.separator+options_intern.folder_name_tgen_preprocessing_binary_trees_sorted);
        f_TGEN_preprocess_binary_trees_sorted.mkdir();

        logger.logLine("[TGENE] Preprocessing GTF.");
        //restructure GTF for TGEN -> only transcript data positions are allowed

        String[] gtf_name_split_dir= options_intern.tepic_gene_annot.split(File.separator);
        String[] gtf_name_split = gtf_name_split_dir[gtf_name_split_dir.length-1].split("\\.");
        String gtf_name = "";
        for(int i = 0; i < gtf_name_split.length-1; i++)
        {
            if(i>0)
            {
                gtf_name+=".";
                gtf_name+=gtf_name_split[i];
            }
            else
            {
                gtf_name+=gtf_name_split[i];

            }
        }

        BufferedReader br_gtf_transcripts = new BufferedReader(new FileReader(new File(options_intern.tepic_gene_annot)));
        BufferedWriter bw_gtf_transcripts = new BufferedWriter(new FileWriter(new File(f_TGEN_preprocess_gtf.getAbsolutePath()+File.separator+gtf_name+options_intern.file_suffix_tgen_preprocess_gtf)));

        String line_gtf_transcripts = "";
        while((line_gtf_transcripts= br_gtf_transcripts.readLine())!=null)
        {
            if(line_gtf_transcripts.startsWith("#"))
            {
                continue;
            }
            String[] split = line_gtf_transcripts.split("\t");
            if(split[2].equals("transcript"))
            {
                if(split[0].matches(".*M.*"))
                {
                    String line = options_intern.tgen_mt_writing;
                    for(int i = 1; i < split.length;i++)
                    {
                        line+="\t"+split[i];
                    }
                    bw_gtf_transcripts.write(line);
                    bw_gtf_transcripts.newLine();
                }
                else
                {
                    if(split[0].matches("chr.*"))
                    {
                        bw_gtf_transcripts.write(line_gtf_transcripts.substring(3));
                        bw_gtf_transcripts.newLine();
                    }
                    else
                    {
                        bw_gtf_transcripts.write(line_gtf_transcripts);
                        bw_gtf_transcripts.newLine();
                    }
                }
            }

        }
        bw_gtf_transcripts.close();
        br_gtf_transcripts.close();

        logger.logLine("[TGENE] Preprocessing Binary Trees for position to Gen Mapping.");

        BufferedReader br_gtf_bin_tree = new BufferedReader(new FileReader(new File(options_intern.tepic_gene_annot)));
        BufferedWriter bw_gtf_bin_tree=new BufferedWriter(new FileWriter(new File(f_TGEN_preprocess_binary_trees_unmerged.getAbsolutePath()+File.separator+"test"+".txt")));

        String current_chr="";
        int count = 0;

        String line_gtf_bin_tree = "";
        while((line_gtf_bin_tree= br_gtf_bin_tree.readLine())!=null)
        {
            if(line_gtf_bin_tree.startsWith("#"))
            {
                continue;
            }
            String[] split = line_gtf_bin_tree.split("\t");
            if(split[2].equals("transcript"))
            {
                String[] split_line = split[8].split(";");

                String chr = split[0];

                if(!current_chr.equals(chr))
                {
                    bw_gtf_bin_tree.close();
                    count=0;
                    bw_gtf_bin_tree=new BufferedWriter(new FileWriter(new File(f_TGEN_preprocess_binary_trees_unmerged.getAbsolutePath()+File.separator+chr+".txt")));
                    bw_gtf_bin_tree.write("#\tPOS_START\tPOS_END\tENSG\n");
                    current_chr=chr;
                }

                String ensg = "";
                String start_pos = split[3];
                String end_pos = split[4];

                for(int i = 0; i < split_line.length; i++)
                {
                    if(split_line[i].startsWith("gene_id"))
                    {
                        String[] split_x = split_line[i].split("\"");
                        ensg=split_x[1].substring(0,split_x[1].length()).split("\\.")[0];
                        break;
                    }
                }

                StringBuilder sb = new StringBuilder();
                sb.append(count);
                sb.append("\t");
                sb.append(start_pos);
                sb.append("\t");
                sb.append(end_pos);
                sb.append("\t");
                sb.append(ensg);

                bw_gtf_bin_tree.write(sb.toString());
                bw_gtf_bin_tree.newLine();

                count++;
            }

        }
        bw_gtf_bin_tree.close();
        br_gtf_bin_tree.close();

        //binary tree merge overlappings

        ArrayList<ENSG_ranges_binary_trees> unmerged_intervalls = new ArrayList<>();

        for(File fileDir: f_TGEN_preprocess_binary_trees_unmerged.listFiles())
        {
            if(!fileDir.isDirectory()&&!fileDir.getName().equals("test.txt"))
            {
                BufferedReader br = new BufferedReader(new FileReader(fileDir));
                String line = br.readLine();
                while((line=br.readLine())!=null)
                {
                    String[] split = line.split("\t");

                    ENSG_ranges_binary_trees iu = new ENSG_ranges_binary_trees();
                    iu.number=Integer.parseInt(split[0]);
                    iu.left_border=Integer.parseInt(split[1]);
                    iu.right_border=Integer.parseInt(split[2]);
                    iu.ensgs.add(split[3]);

                    unmerged_intervalls.add(iu);
                }
                br.close();

                Collections.sort(unmerged_intervalls);

                ArrayList<ENSG_ranges_binary_trees> unmerged_intervalls_temp_list = new ArrayList<>();
                //merge intervalls with same ENSG list!
                int count_same_ENSG = 0;
                boolean last_one_merged=false;
                for(int i = 1; i < unmerged_intervalls.size();i++)
                {
                    ENSG_ranges_binary_trees iu_before = unmerged_intervalls.get(i-1);
                    ENSG_ranges_binary_trees iu_now = unmerged_intervalls.get(i);

                    HashSet<String> ensgs_before = iu_before.ensgs;
                    HashSet<String> ensgs_now = iu_now.ensgs;

                    boolean all_in = true;

                    for(String s: ensgs_before)
                    {
                        if(!ensgs_now.contains(s))
                        {
                            all_in=false;
                        }
                    }
                    for(String s: ensgs_now)
                    {
                        if(!ensgs_before.contains(s))
                        {
                            all_in=false;
                        }
                    }

                    if(all_in)
                    {
                        //now merge!

                        //look for further intervalls with same ensg
                        boolean found_furthers = true;
                        int counter_further = 0;
                        int temp_i = i;
                        while(found_furthers && temp_i < unmerged_intervalls.size())
                        {
                            found_furthers=false;

                            HashSet<String> ensgs_further = unmerged_intervalls.get(temp_i).ensgs;


                            boolean all_in_further = true;

                            for(String s: ensgs_further)
                            {
                                if(!ensgs_before.contains(s))
                                {
                                    all_in_further=false;
                                }
                            }

                            for(String s: ensgs_before)
                            {
                                if(!ensgs_further.contains(s))
                                {
                                    all_in_further=false;
                                }
                            }

                            if(all_in_further)
                            {
                                found_furthers=true;
                                iu_now=unmerged_intervalls.get(temp_i);
                            }
                            else
                            {
                                found_furthers=false;
                                break;
                            }

                            if(all_in_further&&temp_i==unmerged_intervalls.size()-1)
                            {
                                last_one_merged=true;
                            }


                            counter_further++;
                            temp_i++;
                        }
                        i+=counter_further+1;


                        ENSG_ranges_binary_trees iu = new ENSG_ranges_binary_trees();
                        iu.number=count_same_ENSG;
                        iu.ensgs = new HashSet<>(ensgs_before);
                        iu.left_border=iu_before.left_border;
                        iu.right_border=iu_now.right_border;

                        unmerged_intervalls_temp_list.add(iu);

                        count_same_ENSG++;
                    }
                    else
                    {
                        iu_before.number=count_same_ENSG;
                        unmerged_intervalls_temp_list.add(iu_before);
                        count_same_ENSG++;
                    }
                }

                if(!last_one_merged)
                {
                    ENSG_ranges_binary_trees iu = unmerged_intervalls.get(unmerged_intervalls.size()-1);
                    iu.number=unmerged_intervalls_temp_list.size();
                    unmerged_intervalls_temp_list.add(iu);
                }

                unmerged_intervalls.clear();
                unmerged_intervalls = new ArrayList<>(unmerged_intervalls_temp_list);
                unmerged_intervalls_temp_list.clear();

                //re calculate borders

                unmerged_intervalls.get(0).left_border=0;

                for(int i = 1; i < unmerged_intervalls.size(); i++)
                {
                    ENSG_ranges_binary_trees iu_before = unmerged_intervalls.get(i-1);
                    ENSG_ranges_binary_trees iu_current = unmerged_intervalls.get(i);

                    iu_current.left_border=iu_before.right_border+1;
                }

                Collections.sort(unmerged_intervalls);

                for(int i = 0; i < unmerged_intervalls.size();i++)
                {
                    unmerged_intervalls.get(i).number = i;
                }

                BufferedWriter bw = new BufferedWriter(new FileWriter(new File(f_TGEN_preprocess_binary_trees_merged+File.separator+fileDir.getName())));
                bw.write("#\tPOS_START\tPOS_END\tENSG\n");
                for(int i = 1; i < unmerged_intervalls.size(); i++)
                {
                    ENSG_ranges_binary_trees iu = unmerged_intervalls.get(i);
                    StringBuilder sb = new StringBuilder();
                    sb.append(iu.number);
                    sb.append("\t");
                    sb.append(iu.left_border);
                    sb.append("\t");
                    sb.append(iu.right_border);
                    sb.append("\t");
                    sb.append(iu.ensgs_to_String());

                    bw.write(sb.toString());
                    bw.newLine();

                }
                bw.close();
            }
        }



        //sort for binary tree implementation

        logger.logLine("[TGENE] Preparing chromosomes for binary search.");

        for(File fileDir: f_TGEN_preprocess_binary_trees_merged.listFiles())
        {
            if(!fileDir.isDirectory())
            {

                ArrayList<ENSG_ranges_binary_trees> regions = new ArrayList<>();

                BufferedReader br = new BufferedReader(new FileReader(fileDir));
                String header = br.readLine();
                String line = "";

                while((line= br.readLine())!=null)
                {
                    String[] split = line.split("\t");
                    ENSG_ranges_binary_trees iu = new ENSG_ranges_binary_trees();
                    iu.number = Integer.parseInt(split[0]);
                    iu.left_border = Integer.parseInt(split[1]);
                    iu.right_border = Integer.parseInt(split[2]);
                    iu.ensgs.addAll(Arrays.asList(split[3].split(";")));

                    regions.add(iu);
                }
                br.close();

                ArrayList<ENSG_ranges_binary_trees> newly_ordered = new ArrayList<>();

                ENSG_ranges_binary_trees median = regions.get(regions.size()/2);
                newly_ordered.add(median);

                ArrayList<ENSG_ranges_binary_trees> region_left=new ArrayList<>();
                for(int i = 0; i < regions.size()/2; i++)
                {
                    region_left.add(regions.get(i));
                }

                ArrayList<ENSG_ranges_binary_trees> region_right=new ArrayList<>();
                for(int i = regions.size()/2+1; i < regions.size(); i++)
                {
                    region_right.add(regions.get(i));
                }

                newly_ordered = recursive_split(region_left,newly_ordered);
                newly_ordered = recursive_split(region_right,newly_ordered);

                BufferedWriter bw = new BufferedWriter(new FileWriter(f_TGEN_preprocess_binary_trees_sorted.getAbsolutePath()+File.separator+fileDir.getName()));
                bw.write(header);
                bw.newLine();

                for(int i = 0; i < newly_ordered.size();i++)
                {
                    bw.write(newly_ordered.get(i).toString_binary());
                    bw.newLine();
                }

                bw.close();
            }
        }


        logger.logLine("[TGENE] Finished preprocessing.");
    }

    /**
     * postprocesses the TEPIC output (checks for TPM filter and copies files into a structure where preprocessing of DYNAMITE can happen
     */
    public void postprocess_tepic_output() throws Exception
    {
        logger.logLine("Start postprocessing of TEPIC output");

        HashMap<String,HashMap<String,HashSet<String>>> groups_to_compare = checkGroupsTEPIC();

        String suffix="";

        if(options_intern.tepic_tpm_cutoff>0)
        {
            suffix="_Gene_View_Filtered_TPM.txt";
        }
        else
        {
            suffix="_Gene_View_Filtered.txt";
        }

        File folder_postprocessing = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tepic_postprocessing);
        folder_postprocessing.mkdir();

        File folder_pp_input = new File(folder_postprocessing.getAbsolutePath()+File.separator+options_intern.folder_name_tepic_postprocessing_input);
        folder_pp_input.mkdir();

        File folder_pp_output = new File(folder_postprocessing.getAbsolutePath()+File.separator+options_intern.folder_name_tepic_postprocessing_output);
        folder_pp_output.mkdir();

        HashSet<String> available_hms = new HashSet<>();

        HashSet<String> check_tfs = new HashSet<>();
        String seperator= "";

        //generate input structure and copy files
        for(String s : groups_to_compare.keySet())
        {
            File f_input_tp_folders = new File(folder_pp_input.getAbsolutePath()+File.separator+s);
            f_input_tp_folders.mkdir();
            HashMap<String,HashSet<String>> hm = groups_to_compare.get(s);
            for(String ss: hm.keySet())
            {
                available_hms.add(ss);
                File f_input_hm_folders = new File(f_input_tp_folders.getAbsolutePath()+File.separator+ss);
                f_input_hm_folders.mkdir();

                //move coresponding sample outputs to these folders
                File f_input_samples = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tepic_output_raw+File.separator+s+File.separator+ss);

                for(File fileDir : f_input_samples.listFiles())
                {
                    if (fileDir.isDirectory()) {
                        for (File fileDir2 : fileDir.listFiles())
                        {
                            if (!fileDir2.isDirectory()) {
                                String name = fileDir2.getName();
                                if (name.matches(".*" + suffix))
                                {
                                    //check tfs
                                    BufferedReader br = new BufferedReader(new FileReader(fileDir2));
                                    String line = br.readLine();
                                    String[] split = line.split("\t");
                                    for(String st : split)
                                    {
                                        if(st.matches(".*::.*") || st.matches(".*[.][.].*"))
                                        {
                                            check_tfs.add(st);
                                            if(st.matches(".*::.*"))
                                            {
                                                seperator="::";
                                            }
                                            else
                                            {
                                                seperator="[.][.]";
                                            }
                                        }
                                    }

                                    br.close();

                                    //COPY!!
                                    String command = "cp -u " + fileDir2.getAbsolutePath() + " " + folder_pp_input.getAbsolutePath() + File.separator + s + File.separator + ss;
                                    Process child = Runtime.getRuntime().exec(command);
                                    logger.logLine("[TEPIC] Copy files: " + command);
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
                        }
                    }
                }
            }
        }

        File f_out_tfs = new File(folder_postprocessing.getAbsolutePath()+File.separator+options_intern.folder_name_tepic_postprocessing_tfs);
        f_out_tfs.mkdir();
        BufferedWriter bw_check_tfs = new BufferedWriter(new FileWriter(f_out_tfs.getAbsolutePath()+File.separator+options_intern.file_suffix_tepic_postprocessing_tfs_tfs));
        for(String key_tfs : check_tfs)
        {
            String[] split = key_tfs.split(seperator);

            StringBuilder sb = new StringBuilder();

            for(int i = 0; i < split.length;i++)
            {
                if(i<1)
                {
                    sb.append(split[i].toUpperCase());
                }
                else
                {
                    sb.append("..");
                    sb.append(split[i].toUpperCase());
                }
            }

            for(int i = 0; i < split.length;i++)
            {
                sb.append("\t");
                sb.append(split[i].toUpperCase());
            }
            bw_check_tfs.write(sb.toString());
            bw_check_tfs.newLine();
        }
        bw_check_tfs.close();

        //generate output structure
        //HMs
        HashMap<String,File> pp_output_hms_files = new HashMap<>();
        for(String s: available_hms)
        {
            File f = new File(folder_pp_output.getAbsolutePath()+File.separator+s);
            f.mkdir();
            pp_output_hms_files.put(s,f);
        }

        //identify groups to compute mean ratios
        //group1_vs_group2
        HashMap<String,File> pp_output_clashedGroups = new HashMap<>();
        HashSet<String> already_checked_groups = new HashSet<>();
        for(String s: groups_to_compare.keySet())
        {
            for(String ss: groups_to_compare.keySet())
            {
                if(s.equals(ss))
                {
                    continue;
                }
                String key1 = s+"_"+ss;
                String key2 = ss+"_"+s;
                if(already_checked_groups.contains(key1) || already_checked_groups.contains(key2))
                {
                    continue;
                }
                already_checked_groups.add(key1);

                HashMap<String,HashSet<String>> group1_hms = groups_to_compare.get(s);
                HashMap<String,HashSet<String>> group2_hms = groups_to_compare.get(ss);

                for(String k: group1_hms.keySet())
                {
                    if(group2_hms.containsKey(k))
                    {
                        File f = new File(pp_output_hms_files.get(k).getAbsolutePath()+File.separator+key1);
                        f.mkdir();
                    }
                }
            }
        }

        //prepare command line for all (computeMeanRatioTFAffinities.py
        String command_base = "python " + options_intern.path_to_COM2POSE+File.separator+options_intern.directory_for_tepic_DYNAMITE+File.separator+"computeMeanRatioTFAffinities.py";

        for(File fileDir : folder_pp_output.listFiles())
        {
            if(fileDir.isDirectory())
            {
                for(File fileDir2 : fileDir.listFiles())
                {
                    if(fileDir2.isDirectory())
                    {
                        String current_HM = fileDir.getName();
                        String current_group_clash = fileDir2.getName();
                        String[] group_clash_split = current_group_clash.split("_");

                        String group1 = group_clash_split[0];
                        String group2 = group_clash_split[1];

                        String group1_input_dir ="";
                        String group2_input_dir="";
                        //if TPM filter was used we need to postprocess the TPM files. otherwise it will not work!
                        if(options_intern.tepic_tpm_cutoff>0)
                        {
                            logger.logLine("[TEPIC] TPM filter > 0, start postprocessing of TPM filtered scores");
                            File output_post_group1 = new File(folder_pp_output+File.separator+current_HM+File.separator+current_group_clash+File.separator+group1);
                            output_post_group1.mkdir();
                            File output_post_group2 = new File(folder_pp_output+File.separator+current_HM+File.separator+current_group_clash+File.separator+group2);
                            output_post_group2.mkdir();

                            //build intersect and write new files with filter

                            File folder_group1 = new File(folder_pp_input.getAbsolutePath()+File.separator+group1+File.separator+current_HM);
                            File folder_group2 = new File(folder_pp_input.getAbsolutePath()+File.separator+group2+File.separator+current_HM);

                            File[] samples_group1 = folder_group1.listFiles();
                            File[] samples_group2 = folder_group2.listFiles();

                            ArrayList<Boolean> write_group1 = new ArrayList<>();
                            ArrayList<Boolean> write_group2 = new ArrayList<>();

                            ArrayList<String> header_group1;
                            ArrayList<String> header_group2;
                            HashSet<String> header_group1_set;
                            HashSet<String> header_group2_set;

                            BufferedReader br_group1_check = new BufferedReader(new FileReader(samples_group1[0]));
                            String line_group1_check =br_group1_check.readLine();
                            header_group1=new ArrayList<>(Arrays.asList(line_group1_check.split("\t")));
                            header_group1_set=new HashSet<>(Arrays.asList(line_group1_check.split("\t")));
                            br_group1_check.close();

                            BufferedReader br_group2_check = new BufferedReader(new FileReader(samples_group2[0]));
                            String line_group2_check =br_group2_check.readLine();
                            header_group2=new ArrayList<>(Arrays.asList(line_group2_check.split("\t")));
                            header_group2_set=new HashSet<>(Arrays.asList(line_group2_check.split("\t")));
                            br_group2_check.close();

                            for(int i = 0; i < header_group1.size();i++)
                            {
                                if(header_group2_set.contains(header_group1.get(i)))
                                {
                                    write_group1.add(true);
                                }
                                else
                                {
                                    write_group1.add(false);
                                }
                            }
                            for(int i = 0; i < header_group2.size(); i++)
                            {
                                if(header_group1_set.contains(header_group2.get(i)))
                                {
                                    write_group2.add(true);
                                }
                                else
                                {
                                    write_group2.add(false);
                                }
                            }

                            for(File f:folder_group1.listFiles())
                            {
                                BufferedReader br_group1 = new BufferedReader(new FileReader(f));
                                BufferedWriter bw_group1 = new BufferedWriter(new FileWriter(new File(output_post_group1.getAbsolutePath()+File.separator+f.getName())));

                                String line_group1 = "";
                                while((line_group1=br_group1.readLine())!=null)
                                {
                                    String[] split_line_group1 = line_group1.split("\t");
                                    StringBuilder sb = new StringBuilder();

                                    for(int i = 0; i < split_line_group1.length;i++)
                                    {
                                        if(write_group1.get(i))
                                        {
                                            if(i>0)
                                            {
                                                sb.append("\t");
                                                sb.append(split_line_group1[i]);
                                            }
                                            else
                                            {
                                                sb.append(split_line_group1[i]);
                                            }
                                        }
                                    }
                                    bw_group1.write(sb.toString());
                                    bw_group1.newLine();
                                }
                                bw_group1.close();
                                br_group1.close();
                            }


                            for(File f:folder_group2.listFiles())
                            {
                                BufferedReader br_group2 = new BufferedReader(new FileReader(f));
                                BufferedWriter bw_group2 = new BufferedWriter(new FileWriter(new File(output_post_group2.getAbsolutePath()+File.separator+f.getName())));

                                String line_group2 = "";
                                while((line_group2=br_group2.readLine())!=null)
                                {
                                    String[] split_line_group2 = line_group2.split("\t");
                                    StringBuilder sb = new StringBuilder();

                                    for(int i = 0; i < split_line_group2.length;i++)
                                    {
                                        if(write_group2.get(i))
                                        {
                                            if(i>0)
                                            {
                                                sb.append("\t");
                                                sb.append(split_line_group2[i]);
                                            }
                                            else
                                            {
                                                sb.append(split_line_group2[i]);
                                            }
                                        }
                                    }
                                    bw_group2.write(sb.toString());
                                    bw_group2.newLine();
                                }
                                bw_group2.close();
                                br_group2.close();
                            }
                            logger.logLine("[TEPIC] TPM filter > 0, end postprocessing of TPM filtered scores");

                            //set to postprocessed TMP filtered data
                            group1_input_dir=output_post_group1.getAbsolutePath();
                            group2_input_dir=output_post_group2.getAbsolutePath();
                        }
                        else
                        {
                            group1_input_dir = folder_pp_input.getAbsolutePath()+File.separator+group1+File.separator+current_HM;
                            group2_input_dir = folder_pp_input.getAbsolutePath()+File.separator+group2+File.separator+current_HM;
                        }


                        File output_mean_affinities = new File(folder_pp_output+File.separator+current_HM+File.separator+current_group_clash+File.separator+options_intern.folder_name_tepic_postprocessing_output_mean_affinities);
                        output_mean_affinities.mkdir();
                        File output_ratios = new File(folder_pp_output+File.separator+current_HM+File.separator+current_group_clash+File.separator+options_intern.folder_name_tepic_postprocessing_output_ratios);
                        output_ratios.mkdir();

                        File output_mean_affinities_group1 = new File(output_mean_affinities.getAbsolutePath()+File.separator+options_intern.file_suffix_tepic_postprocessing_output_mean_affinities+group1+".txt");
                        File output_mean_affinities_group2 = new File(output_mean_affinities.getAbsolutePath()+File.separator+options_intern.file_suffix_tepic_postprocessing_output_mean_affinities+group2+".txt");
                        File output_ratios_group1_group2 = new File(output_ratios.getAbsolutePath()+File.separator+options_intern.file_suffix_tepic_postprocessing_output_ratios+group1+"_"+group2+".txt");

                        HashSet<File> files_to_create = new HashSet<>();
                        files_to_create.add(output_mean_affinities_group1);
                        files_to_create.add(output_mean_affinities_group2);
                        files_to_create.add(output_ratios_group1_group2);

                        for(File f: files_to_create)
                        {
                            BufferedWriter bw = new BufferedWriter(new FileWriter(f));
                            bw.write("");
                            bw.close();
                        }


                        String command_edited = new String(command_base);

                        command_edited+= " " + group1_input_dir+File.separator;
                        command_edited+= " " + group2_input_dir+File.separator;
                        command_edited+= " " + output_mean_affinities_group1.getAbsolutePath();
                        command_edited+= " " + output_mean_affinities_group2.getAbsolutePath();
                        command_edited+= " " + output_ratios_group1_group2.getAbsolutePath();

                        String command_tail = "";

                        if(options_intern.tepic_original_decay)
                        {
                            command_tail += " True";
                        }
                        else
                        {
                            command_tail += " False";
                        }
                        if(options_intern.tepic_not_generated)
                        {
                            command_tail += " True";
                        }
                        else
                        {
                            command_tail+= " False";
                        }
                        if(options_intern.tepic_tpm_cutoff>0)
                        {
                            command_tail+= " True";
                        }
                        else
                        {
                            command_tail+= " False";
                        }

                        command_edited += command_tail;

                        logger.logLine("[TEPIC] execute computeMeanRatioTFAffinities.py with command line: " + command_edited);
                        Process child = Runtime.getRuntime().exec(command_edited);
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
            }
        }
        logger.logLine("Finished postprocessing of TEPIC output");
    }

    /**
     * creates the TEPIC command lines and runs all samples of all histone modification and all timepoints
     * @throws Exception
     */
    public void run_tepic() throws Exception {
        logger.logLine("Start TEPIC.sh");

        String command = "bash";
        String tepic_path = " " + options_intern.path_to_COM2POSE+File.separator+options_intern.directory_for_tepic_scripts_code_tepic_sh;
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

        File output_TEPIC = new File(options_intern.com2pose_working_directory+File.separator+ options_intern.folder_name_tepic_output_raw);
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
                            String n_dir = options_intern.com2pose_working_directory+File.separator+ options_intern.folder_name_deseq2_preprocessing+File.separator+options_intern.folder_name_deseq2_preprocessing_single+File.separator+dirGroup.getName()+options_intern.file_suffix_deseq2_preprocessing_meanCounts;
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

        File folder = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_deseq2_R_scripts);

        for (File dir:folder.listFiles())
        {
            if(!dir.isDirectory())
            {
                String command = "Rscript " + dir.getAbsolutePath();
                Process child = Runtime.getRuntime().exec(command);
                logger.logLine("[DESEQ2] Running script " + dir.getName()+": " + command);
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


        File folder_results = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_deseq2_output_raw);
        File output_file = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_deseq2_output);
        output_file.mkdir();

        for(File res:folder_results.listFiles())
        {
            if(!res.isDirectory())
            {
                String name_res = res.getName();
                String[] split_name_res = name_res.split("_");
                String name_out = split_name_res[0]+"_"+split_name_res[1]+options_intern.file_suffix_deseq2_output_DYNAMITE;

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



        logger.logLine("Finished postprocessing DESeq2 data for input to DYNAMITE");

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

        File output_intermediate_steps = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_deseq2_preprocessing);
        if(output_intermediate_steps.exists())
        {
            //logger.logLine("Working directory was already used - please us another one or empty this one completely");
            //TODO: after debugging use system exit !!!
            //System.exit(1);
        }
        output_intermediate_steps.mkdir();
        File output_inter_steps_combined=new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_deseq2_preprocessing+File.separator+options_intern.folder_name_deseq2_preprocessing_combined);
        output_inter_steps_combined.mkdir();
        File output_inter_steps_single=new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_deseq2_preprocessing+File.separator+options_intern.folder_name_deseq2_preprocessing_single);
        output_inter_steps_single.mkdir();
        File output_inter_steps_symbols_ensg_mean_counts = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_deseq2_preprocessing+File.separator+options_intern.folder_name_deseq2_preprocessing_gene_symbols);
        output_inter_steps_symbols_ensg_mean_counts.mkdir();

        File folder = new File(options_intern.deseq2_input_directory);

        HashMap<String,String> ensg_symbol = new HashMap<>();
        BufferedReader br_ensg_symbol = new BufferedReader(new FileReader(new File(options_intern.tepic_ensg_symbol)));
        String line_ensg_symbol = br_ensg_symbol.readLine();
        while((line_ensg_symbol=br_ensg_symbol.readLine())!=null)
        {
            String[] split = line_ensg_symbol.split("\t");
            if(split.length>1)
            {
                ensg_symbol.put(split[0],split[1]);
            }
        }

        br_ensg_symbol.close();

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

                BufferedWriter bw_means = new BufferedWriter(new FileWriter(new File(output_inter_steps_single.getAbsolutePath()+File.separator+group.getName()+options_intern.file_suffix_deseq2_preprocessing_meanCounts)));
                bw_means.write(group.getName()+"_MEANS");
                bw_means.newLine();
                for(int i = 0; i < count_ensg_lines; i++)
                {
                    int mean_count = mean_line_counts.get(i)/count_samples;
                    bw_means.write(""+mean_count);
                    bw_means.newLine();

                }
                bw_means.close();

                BufferedReader br_ensg = new BufferedReader(new FileReader(new File(options_intern.deseq2_input_gene_id)));
                BufferedWriter bw_symbol = new BufferedWriter(new FileWriter(output_inter_steps_symbols_ensg_mean_counts.getAbsolutePath()+File.separator+group.getName()+".csv"));
                bw_symbol.write("SYMBOL\tENSG\tMEAN_COUNT");
                bw_symbol.newLine();

                String line = br_ensg.readLine();

                int i = 0;
                while((line=br_ensg.readLine())!=null)
                {
                    StringBuilder sb = new StringBuilder();

                    int mean_count = mean_line_counts.get(i)/count_samples;
                    String gene_symbol_name = "NO_SYMBOL";

                    if(ensg_symbol.containsKey(line))
                    {
                        gene_symbol_name = ensg_symbol.get(line);
                    }

                    sb.append(gene_symbol_name);
                    sb.append("\t");
                    sb.append(line);
                    sb.append("\t");
                    sb.append(mean_count);


                    bw_symbol.write(sb.toString());
                    bw_symbol.newLine();
                    i++;
                }
                br_ensg.close();
                bw_symbol.close();



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
        File r_scripts = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_deseq2_R_scripts);
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

                File output_deseq2 = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_deseq2_output_raw);
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

    public void get_ensg_symbol_mapping() throws Exception {

        logger.logLine("[DESEQ2] Start mapping ENSG to GENE SYMBOLS.");

        File script_output_dir = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_deseq2_preprocessing);
        script_output_dir.mkdir();

        File results = new File(script_output_dir.getAbsolutePath()+File.separator+options_intern.file_suffix_deseq2_mapping);

        StringBuilder sb = new StringBuilder();
        sb.append("if (!requireNamespace(\"BiocManager\", quietly = TRUE))\n" +
                "  install.packages(\"BiocManager\")\n" +
                "\n" +
                "if (!requireNamespace(\"biomaRt\", quietly = TRUE))\n" +
                "  install.packages(\"biomaRt\")\n" +
                "\n" +
                "library('biomaRt')\n");

        sb.append("df <- read.csv('"+options_intern.deseq2_input_gene_id+"')\n");

        sb.append("mart <- useDataset(\""+options_intern.deseq2_biomart_dataset_species+"\", useMart(\"ensembl\"))");
        sb.append("\n");

        sb.append("df$id <- NA\n" +
                "G_list <- getBM(filters= \"ensembl_gene_id\", attributes= c(\"ensembl_gene_id\",\"mgi_symbol\"),values=df$Geneid,mart= mart)\n" +
                "write.table(G_list,\""+results.getAbsolutePath()+"\", row.names = FALSE, quote = F, sep=\"\\t\")\n");


        File script_output = new File(script_output_dir.getAbsolutePath()+File.separator+"ENSG_SYMBOL_MAP.R");

        BufferedWriter bw = new BufferedWriter(new FileWriter(script_output));
        bw.write(sb.toString());
        bw.close();

        String command = "Rscript " + script_output.getAbsolutePath();
        Process child = Runtime.getRuntime().exec(command);
        logger.logLine("[DESEQ2] Running script " + script_output.getName()+": " + command);
        int code = child.waitFor();
        switch (code) {
            case 0:
                break;
            case 1:
                String message = child.getErrorStream().toString();
                throw new Exception(message);
        }

        options_intern.tepic_ensg_symbol = results.getAbsolutePath();


        logger.logLine("[DESEQ2] Finished mapping ENSG to GENE SYMBOLS.");

    }

    /**
     * Filter TEPIC input files for blacklisted regions
     */
    public void filter_blacklist() throws IOException {
        File folder_input = new File(options_intern.tepic_input_directory);
        File output_folder = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_blacklisted_regions);
        File output_folder_new_input = new File(output_folder.getAbsolutePath()+File.separator+options_intern.folder_name_blacklisted_regions_new_input);
        output_folder_new_input.mkdir();

        //set new folder directory for tepic input and save old one
        options_intern.tepic_input_prev=options_intern.tepic_input_directory;
        options_intern.tepic_input_directory = output_folder_new_input.getAbsolutePath();

        logger.logLine("[BLACKLIST] Create chromosome binary trees.");
        //CREATE BINARY TREES
        HashMap<String,BL_binary_tree> chr_tree = new HashMap<>();

        File input_folder_chr = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_blacklisted_regions+File.separator+options_intern.folder_name_blacklisted_regions_preprocessing+File.separator+options_intern.folder_name_blacklisted_regions_preprocessing_sorted);
        for(File fileChr : input_folder_chr.listFiles())
        {
            if(fileChr.isFile())
            {
                String name= fileChr.getName().split("\\.")[0];

                ArrayList<BL_ranges_binary_tree> region = new ArrayList<>();

                BufferedReader br_chr = new BufferedReader(new FileReader(fileChr));
                String line_chr = br_chr.readLine();
                while((line_chr=br_chr.readLine())!=null)
                {
                    String[] split = line_chr.split("\t");

                    BL_ranges_binary_tree iu = new BL_ranges_binary_tree();
                    iu.number=Integer.parseInt(split[0]);
                    iu.chr = split[1];
                    iu.left_border = Integer.parseInt(split[2]);
                    iu.right_border = Integer.parseInt(split[3]);
                    iu.signal = split[4];

                    region.add(iu);
                }
                br_chr.close();

                BL_binary_tree_node root = new BL_binary_tree_node(region.get(0),region.get(0).number);
                BL_binary_tree tree = new BL_binary_tree(root);

                for(int i = 1; i < region.size(); i++)
                {
                    tree.add(region.get(i).number,region.get(i));
                }

                chr_tree.put(name,tree);

            }
        }

        logger.logLine("[BLACKLIST] Filter input files for blacklisted regions.");
        //now filter all files
        int count_matched = 0;

        for(File fileDirTP : folder_input.listFiles())
        {
            if(fileDirTP.isDirectory())
            {
                File output_folder_new_input_TP = new File(output_folder_new_input.getAbsolutePath()+File.separator+fileDirTP.getName());
                output_folder_new_input_TP.mkdir();

                for(File fileDirTP_HM : fileDirTP.listFiles())
                {
                    if(fileDirTP_HM.isDirectory())
                    {
                        File output_folder_new_input_TP_HM = new File(output_folder_new_input_TP.getAbsolutePath()+File.separator+fileDirTP_HM.getName());
                        output_folder_new_input_TP_HM.mkdir();

                        for(File filrDirTP_HM_sample : fileDirTP_HM.listFiles())
                        {
                            if(filrDirTP_HM_sample.isFile())
                            {
                                logger.logLine("[BLACKLIST] Filter: " + fileDirTP.getName() + ": " + fileDirTP_HM.getName() + " - " + filrDirTP_HM_sample.getName());
                                //FILTER HERE WITH BINARY TREE

                                BufferedWriter bw = new BufferedWriter(new FileWriter(new File(output_folder_new_input_TP_HM.getAbsolutePath()+File.separator+filrDirTP_HM_sample.getName())));
                                BufferedReader br = new BufferedReader(new FileReader(filrDirTP_HM_sample));
                                String line ="";
                                int count_line = 0;
                                while((line=br.readLine())!=null)
                                {
                                    String[] split = line.split("\t");

                                    String chr = split[0];
                                    if(!chr.matches(".*chr.*"))
                                    {
                                        chr="chr"+chr;
                                    }

                                    BL_ranges_binary_tree iu = new BL_ranges_binary_tree();
                                    iu.chr = chr;
                                    iu.left_border = Integer.parseInt(split[1]);
                                    iu.right_border = Integer.parseInt(split[2]);

                                    if(!chr_tree.containsKey(chr))
                                    {
                                        count_line++;
                                        continue;
                                    }

                                    BL_binary_tree tree = chr_tree.get(chr);
                                    if(tree.containsNode(iu) == null)
                                    {
                                        bw.write(line);
                                        bw.newLine();
                                    }
                                    else
                                    {
                                        count_matched++;
                                    }
                                    count_line++;
                                }
                                br.close();
                                bw.close();
                            }
                        }
                    }
                }
            }
        }
        logger.logLine("[BLACKLIST] Finished filtering for blacklisted regions.");
    }

    /**
     * preprocesses the blacklist file for binary search
     */
    public void preprocess_blacklist() throws IOException {

        logger.logLine("[BLACKLIST] start preprocessing blacklist");

        File f_blacklist = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_blacklisted_regions);
        f_blacklist.mkdir();
        File f_blacklist_pre = new File(f_blacklist.getAbsolutePath()+File.separator+options_intern.folder_name_blacklisted_regions_preprocessing);
        f_blacklist_pre.mkdir();
        File f_blacklist_pre_chr = new File(f_blacklist_pre.getAbsolutePath()+File.separator+options_intern.folder_name_blacklisted_regions_preprocessing_perChr);
        f_blacklist_pre_chr.mkdir();

        BufferedReader br_chr = new BufferedReader(new FileReader(new File(options_intern.black_list_dir)));
        ArrayList<BL_ranges_binary_tree> chr_ius= new ArrayList<>();

        String line_chr = "";
        String current_chr = "";
        BufferedWriter bw_chr = new BufferedWriter(new FileWriter(new File(f_blacklist_pre_chr.getAbsolutePath()+File.separator+"test.txt")));
        while((line_chr=br_chr.readLine())!=null)
        {
            String[] split = line_chr.split("\t");

            String chr = split[0];
            if(!chr.equals(current_chr))
            {
                Collections.sort(chr_ius);

                int i = 0;

                for(BL_ranges_binary_tree iu : chr_ius)
                {
                    iu.number=i;
                    bw_chr.write(iu.toString());
                    bw_chr.newLine();
                    i++;
                }

                chr_ius.clear();

                if(!chr.matches(".*chr.*"))
                {
                    chr="chr"+chr;
                }

                bw_chr.close();
                bw_chr = new BufferedWriter(new FileWriter(new File(f_blacklist_pre_chr.getAbsolutePath()+File.separator+chr+".txt")));
                bw_chr.write("#\tCHR\tLEFT_BORDER\tRIGHT_BORDER\tSIGNAL");
                bw_chr.newLine();
                current_chr=split[0];
            }

            BL_ranges_binary_tree iu = new BL_ranges_binary_tree();
            iu.chr=chr;
            iu.left_border = Integer.parseInt(split[1]);
            iu.right_border = Integer.parseInt(split[2]);
            iu.signal = split[3].replace(' ', '_').toUpperCase();

            if(options_intern.black_list_signals.contains(iu.signal))
            {
                chr_ius.add(iu);
            }
        }

        Collections.sort(chr_ius);

        int i = 0;

        for(BL_ranges_binary_tree iu : chr_ius)
        {
            iu.number=i;
            bw_chr.write(iu.toString());
            bw_chr.newLine();
            i++;
        }

        bw_chr.close();
        br_chr.close();

        logger.logLine("[BLACKLIST] prepare chromosomes for binary tree");

        File f_blacklist_pre_sorted = new File(f_blacklist_pre.getAbsolutePath()+File.separator+options_intern.folder_name_blacklisted_regions_preprocessing_sorted);
        f_blacklist_pre_sorted.mkdir();

        for(File fileDir : f_blacklist_pre_chr.listFiles())
        {
            if(fileDir.isFile() && !fileDir.getName().equals("test.txt"))
            {
                ArrayList<BL_ranges_binary_tree> current_ius = new ArrayList<>();
                String header;

                BufferedReader br_sort = new BufferedReader(new FileReader(fileDir));
                header= br_sort.readLine();
                String line_sort ="";
                while((line_sort=br_sort.readLine())!=null)
                {
                    String[] split = line_sort.split("\t");

                    BL_ranges_binary_tree iu = new BL_ranges_binary_tree();
                    iu.number = Integer.parseInt(split[0]);
                    iu.chr = split[1];
                    iu.left_border = Integer.parseInt(split[2]);
                    iu.right_border = Integer.parseInt(split[3]);
                    iu.signal=split[4];

                    current_ius.add(iu);
                }
                br_sort.close();

                ArrayList<BL_ranges_binary_tree> newly_ordered = new ArrayList<>();

                recursive_split_BL(current_ius,newly_ordered);

                BufferedWriter bw_sort = new BufferedWriter(new FileWriter(f_blacklist_pre_sorted.getAbsolutePath()+File.separator+fileDir.getName()));
                bw_sort.write(header);
                bw_sort.newLine();

                for(int j = 0; j < newly_ordered.size(); j++)
                {
                    bw_sort.write(newly_ordered.get(j).toString());
                    bw_sort.newLine();
                }
                bw_sort.close();
            }

        }


        logger.logLine("[BLACKLIST] finished preprocessing blacklist");


    }

    /**
     * preprocess mix histones, search for same peaks and use either the union or the intersection of all
     */
    public void mix_option() throws IOException {
        File file_root_input = new File(options_intern.tepic_input_directory);
        File root_mix_working_dir = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_mix_option);
        root_mix_working_dir.mkdir();

        File f_sample_mix_preprocess = new File(root_mix_working_dir.getAbsolutePath()+File.separator+options_intern.folder_name_mix_option_sample_mix_preprocessing);
        f_sample_mix_preprocess.mkdir();

        File f_sample_mix_output = new File(root_mix_working_dir.getAbsolutePath()+File.separator+options_intern.folder_name_mix_option_sample_mix);
        f_sample_mix_output.mkdir();

        if(options_intern.mix_level.equals("SAMPLE_LEVEL"))
        {
            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory=f_sample_mix_output.getAbsolutePath();
        }

        logger.logLine("[MIX] Preprocess input data for sample mix - split chromosomes.");

        for(File fileDir:file_root_input.listFiles())
        {
            if(fileDir.isDirectory())
            {
                String timepoint = fileDir.getName();
                File file_output_tp = new File(f_sample_mix_preprocess + File.separator + timepoint);
                file_output_tp.mkdir();

                for(File fileDirHM : fileDir.listFiles())
                {
                    if (fileDirHM.isDirectory())
                    {
                        File file_output_tp_hm = new File(file_output_tp.getAbsolutePath() + File.separator + fileDirHM.getName());
                        file_output_tp_hm.mkdir();

                        for(File fileDirHM_sample : fileDirHM.listFiles())
                        {
                            if(fileDirHM_sample.isFile())
                            {
                                BufferedWriter bw = new BufferedWriter(new FileWriter(new File(file_output_tp_hm.getAbsolutePath()+File.separator+"test.txt")));
                                BufferedReader br = new BufferedReader(new FileReader(fileDirHM_sample));
                                String line = "";
                                String currentChr ="";
                                while((line=br.readLine())!=null)
                                {
                                    String[] split = line.split("\t");
                                    if(!split[0].equals(currentChr))
                                    {
                                        bw.close();
                                        currentChr=split[0];
                                        File f_output_chr = new File(file_output_tp_hm.getAbsolutePath()+File.separator+currentChr);
                                        f_output_chr.mkdir();

                                        bw = new BufferedWriter(new FileWriter(new File(f_output_chr.getAbsolutePath()+File.separator+fileDirHM_sample.getName())));
                                    }
                                    bw.write(line);
                                    bw.newLine();
                                }
                                bw.close();
                            }
                        }
                    }
                }
            }
        }

        logger.logLine("[MIX] Create "+options_intern.mix_option+" of samples");

        for(File fileDir: f_sample_mix_preprocess.listFiles())
        {
            if(fileDir.isDirectory())
            {
                String timepoint = fileDir.getName();
                File f_output_union_samples_tp = new File(f_sample_mix_output.getAbsolutePath()+File.separator+timepoint);
                f_output_union_samples_tp.mkdir();

                for(File fileDirHM : fileDir.listFiles())
                {
                    if(fileDirHM.isDirectory())
                    {
                        String hm = fileDirHM.getName();
                        File f_output_union_samples_tp_hm = new File(f_output_union_samples_tp.getAbsolutePath()+File.separator+hm);
                        f_output_union_samples_tp_hm.mkdir();

                        String file_ending = "";

                        HashMap<String,ArrayList<MIX_Interval>> all_chromosomes = new HashMap();
                        ArrayList<Integer> chr_alpha = new ArrayList<>();
                        ArrayList<String> chr_str = new ArrayList<>();

                        for(File fileDirHM_Chr : fileDirHM.listFiles())
                        {
                            if(fileDirHM_Chr.isDirectory())
                            {
                                String chr = fileDirHM_Chr.getName();

                                ArrayList<MIX_Interval> all_intervals = new ArrayList<>();

                                int sample_number = 0;

                                for(File fileDirHM_Chr_sample : fileDirHM_Chr.listFiles())
                                {
                                    if(fileDirHM_Chr_sample.isFile())
                                    {
                                        String[] f_ending = fileDirHM_Chr_sample.getName().split("\\.");
                                        file_ending=f_ending[f_ending.length-1];
                                        //build array of a union of all peaks
                                        BufferedReader br = new BufferedReader(new FileReader(fileDirHM_Chr_sample));
                                        String line = "";
                                        while((line=br.readLine())!=null)
                                        {
                                            String[] split = line.split("\t");

                                            MIX_Interval_Object mio = new MIX_Interval_Object(split[0],Integer.parseInt(split[1]),Integer.parseInt(split[2]),split[3],Integer.parseInt(split[4]),split[5],Double.parseDouble(split[6]),Double.parseDouble(split[7]),Double.parseDouble(split[8]));
                                            MIX_Interval mi = new MIX_Interval(Integer.parseInt(split[1]),Integer.parseInt(split[2]));
                                            mi.merged_intervals.add(mio);

                                            all_intervals.add(mi);

                                        }
                                        br.close();
                                        sample_number++;

                                        /*
                                        if(options_intern.mix_option.equals("INTERSECTION"))
                                        {
                                            ArrayList<MIX_Interval> intersec_al = new ArrayList<>();

                                            BufferedReader br_intersec = new BufferedReader(new FileReader(fileDirHM_Chr_sample));
                                            String line_intersec = "";
                                            while((line_intersec=br_intersec.readLine())!=null)
                                            {
                                                String[] split = line_intersec.split("\t");

                                                MIX_Interval_Object mio = new MIX_Interval_Object(split[0],Integer.parseInt(split[1]),Integer.parseInt(split[2]),split[3],Integer.parseInt(split[4]),split[5],Double.parseDouble(split[6]),Double.parseDouble(split[7]),Double.parseDouble(split[8]));
                                                MIX_Interval mi = new MIX_Interval(Integer.parseInt(split[1]),Integer.parseInt(split[2]));
                                                mi.merged_intervals.add(mio);

                                                intersec_al.add(mi);

                                            }
                                            br.close();


                                            intersec_vals.put(fileDirHM_Chr_sample.getName(),intersec_al);
                                        }*/

                                    }
                                }

                                //intersections sorting!
                                /*if(options_intern.mix_option.equals("INTERSECTION"))
                                {
                                    for(String k: intersec_vals.keySet())
                                    {
                                        Collections.sort(intersec_vals.get(k));
                                    }
                                }*/

                                //unions sorting
                                Collections.sort(all_intervals);

                                Stack<MIX_Interval> stack_union = mergeIntervals(all_intervals);

                                ArrayList<MIX_Interval> chr_unions = new ArrayList<>();

                                int min_occurence = 0;
                                if(options_intern.mix_occurence_intersection==-1)
                                {
                                    min_occurence=sample_number;
                                }
                                else
                                {
                                    if(sample_number<options_intern.mix_occurence_intersection)
                                    {
                                        min_occurence=sample_number;

                                    }
                                    else
                                    {
                                        min_occurence=options_intern.mix_occurence_intersection;
                                    }
                                }

                                while(!stack_union.isEmpty())
                                {
                                    MIX_Interval t = stack_union.pop();
                                    t.calculate_mean("SAMPLE_LEVEL");

                                    if(options_intern.mix_option.equals("INTERSECTION"))
                                    {
                                        if(t.merged_intervals.size()>=min_occurence)
                                        {
                                            chr_unions.add(t);
                                        }
                                    }
                                    else
                                    {
                                        chr_unions.add(t);
                                    }
                                }

                                Collections.sort(chr_unions);

                                all_chromosomes.put(chr,chr_unions);

                                try {
                                    chr_alpha.add(Integer.parseInt(chr));
                                }
                                catch (Exception e)
                                {
                                    chr_str.add(chr);
                                }
                            }
                        }

                        Collections.sort(chr_alpha);
                        Collections.sort(chr_str);


                        //print Sample_unions => can be used if not HM_Level is used
                        BufferedWriter bw = new BufferedWriter(new FileWriter(f_output_union_samples_tp_hm.getAbsolutePath()+File.separator+timepoint+"_"+hm+"."+file_ending));

                        int peak_counter = 1;

                        for(int chr : chr_alpha)
                        {
                            ArrayList<MIX_Interval> x = all_chromosomes.get(""+chr);
                            for(int i = 0; i < x.size();i++)
                            {
                                bw.write(x.get(i).meanToString(peak_counter));
                                bw.newLine();
                                peak_counter++;
                            }
                        }

                        for(String chr: chr_str)
                        {
                            ArrayList<MIX_Interval> x = all_chromosomes.get(chr);
                            for(int i = 0; i < x.size();i++)
                            {
                                bw.write(x.get(i).meanToString(peak_counter));
                                bw.newLine();
                                peak_counter++;
                            }
                        }
                        bw.close();


                    }
                }
            }
        }

        if(options_intern.mix_level.equals("HM_LEVEL"))
        {

            logger.logLine("[MIX] Preprocess sample unions for HM mix");

            File f_output_preprocessing_hm = new File(root_mix_working_dir.getAbsolutePath()+File.separator+options_intern.folder_name_mix_option_preprocess_hm_mix);
            f_output_preprocessing_hm.mkdir();

            File f_output_hm = new File(root_mix_working_dir.getAbsolutePath()+File.separator+options_intern.folder_name_mix_option_hm_mix);
            f_output_hm.mkdir();

            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory=f_output_hm.getAbsolutePath();

            logger.logLine("[MIX] Identify possible Timepoints with same Histone Modifications");

            HashMap<String,ArrayList<String>> timepoints_histone_modifications = new HashMap<>();
            HashMap<String,ArrayList<String>> deleted_tps = new HashMap<>();
            HashSet<String> available_hms = new HashSet<>();

            //identify possible timepoints
            for(File fileDir:f_sample_mix_output.listFiles())
            {
                if(fileDir.isDirectory())
                {
                    String timepoint = fileDir.getName();
                    ArrayList hm = new ArrayList();

                    for(File fileDirHM : fileDir.listFiles())
                    {
                        if(fileDirHM.isDirectory())
                        {
                            hm.add(fileDirHM.getName());
                            available_hms.add(fileDirHM.getName());
                        }
                    }
                    timepoints_histone_modifications.put(timepoint,hm);
                }
            }


            for(String tp : timepoints_histone_modifications.keySet())
            {
                if(timepoints_histone_modifications.get(tp).size()< available_hms.size())
                {
                    deleted_tps.put(tp,timepoints_histone_modifications.get(tp));
                    timepoints_histone_modifications.remove(tp);
                    continue;
                }

                boolean all_in = true;
                for(String hm : available_hms)
                {
                    if(!timepoints_histone_modifications.get(tp).contains(hm))
                    {
                        all_in=false;
                    }
                }

                if(!all_in)
                {
                    deleted_tps.put(tp,timepoints_histone_modifications.get(tp));
                    timepoints_histone_modifications.remove(tp);
                }
            }

            StringBuilder sb_found_mixing_tps = new StringBuilder();
            sb_found_mixing_tps.append("[MIX] Can perform complete mix for HMs (");
            for(String hm : available_hms)
            {
                sb_found_mixing_tps.append(hm);
                sb_found_mixing_tps.append(" ");
            }
            sb_found_mixing_tps.append(") in timepoints (");
            for(String tp:timepoints_histone_modifications.keySet())
            {
                sb_found_mixing_tps.append(tp);
                sb_found_mixing_tps.append(" ");
            }
            sb_found_mixing_tps.append("). Can perform part-mix or no-mix for timepoints (");
            for(String tp : deleted_tps.keySet())
            {
                sb_found_mixing_tps.append(tp);
                sb_found_mixing_tps.append(" ");
            }
            sb_found_mixing_tps.append(").");
            logger.logLine(sb_found_mixing_tps.toString());


            for(File fileDir : f_sample_mix_output.listFiles())
            {
                if(fileDir.isDirectory())
                {
                    String timepoint = fileDir.getName();
                    File f_output_hm_prepro = new File(f_output_preprocessing_hm.getAbsolutePath()+File.separator+timepoint);
                    f_output_hm_prepro.mkdir();

                    for(File fileDirHM : fileDir.listFiles())
                    {
                        if(fileDirHM.isDirectory())
                        {
                            String hm = fileDirHM.getName();
                            File f_output_hm_prepro_hm = new File(f_output_hm_prepro.getAbsolutePath()+File.separator+"MIX");
                            f_output_hm_prepro_hm.mkdir();

                            for(File fileDirHM_samples : fileDirHM.listFiles())
                            {
                                if(fileDirHM_samples.isFile())
                                {
                                    BufferedWriter bw = new BufferedWriter(new FileWriter(new File(f_output_hm_prepro_hm.getAbsolutePath()+File.separator+"test.txt")));
                                    BufferedReader br = new BufferedReader(new FileReader(fileDirHM_samples));
                                    String line = "";
                                    String currentChr ="";
                                    while((line=br.readLine())!=null)
                                    {
                                        String[] split = line.split("\t");
                                        if(!split[0].equals(currentChr))
                                        {
                                            bw.close();
                                            currentChr=split[0];
                                            File f_output_chr = new File(f_output_hm_prepro_hm.getAbsolutePath()+File.separator+currentChr);
                                            f_output_chr.mkdir();

                                            bw = new BufferedWriter(new FileWriter(new File(f_output_chr.getAbsolutePath()+File.separator+fileDirHM_samples.getName())));
                                        }
                                        bw.write(line);
                                        bw.newLine();
                                    }
                                    bw.close();
                                }
                            }
                        }
                    }
                }
            }

            logger.logLine("[MIX] Create "+options_intern.mix_option+" of HMs");

            for(File fileDir: f_output_preprocessing_hm.listFiles())
            {
                if(fileDir.isDirectory())
                {
                    String timepoint = fileDir.getName();
                    File f_output_union_samples_tp = new File(f_output_hm.getAbsolutePath()+File.separator+timepoint);
                    f_output_union_samples_tp.mkdir();

                    for(File fileDirHM : fileDir.listFiles())
                    {
                        if(fileDirHM.isDirectory())
                        {
                            String hm = fileDirHM.getName();
                            File f_output_union_samples_tp_hm = new File(f_output_union_samples_tp.getAbsolutePath()+File.separator+hm);
                            f_output_union_samples_tp_hm.mkdir();

                            String file_ending = "";

                            HashMap<String,ArrayList<MIX_Interval>> all_chromosomes = new HashMap();
                            ArrayList<Integer> chr_alpha = new ArrayList<>();
                            ArrayList<String> chr_str = new ArrayList<>();

                            for(File fileDirHM_Chr : fileDirHM.listFiles())
                            {
                                if(fileDirHM_Chr.isDirectory())
                                {
                                    String chr = fileDirHM_Chr.getName();

                                    ArrayList<MIX_Interval> all_intervals = new ArrayList<>();

                                    int sample_number = 0;

                                    for(File fileDirHM_Chr_sample : fileDirHM_Chr.listFiles())
                                    {
                                        if(fileDirHM_Chr_sample.isFile())
                                        {
                                            String[] f_ending = fileDirHM_Chr_sample.getName().split("\\.");
                                            file_ending=f_ending[f_ending.length-1];
                                            //build array of a union of all peaks
                                            BufferedReader br = new BufferedReader(new FileReader(fileDirHM_Chr_sample));
                                            String line = "";
                                            while((line=br.readLine())!=null)
                                            {
                                                String[] split = line.split("\t");

                                                MIX_Interval_Object mio = new MIX_Interval_Object(split[0],Integer.parseInt(split[1]),Integer.parseInt(split[2]),split[3],Integer.parseInt(split[4]),split[5],Double.parseDouble(split[6]),Double.parseDouble(split[7]),Double.parseDouble(split[8]));
                                                MIX_Interval mi = new MIX_Interval(Integer.parseInt(split[1]),Integer.parseInt(split[2]));
                                                mi.merged_intervals.add(mio);

                                                all_intervals.add(mi);

                                            }
                                            br.close();
                                            sample_number++;

                                        /*
                                        if(options_intern.mix_option.equals("INTERSECTION"))
                                        {
                                            ArrayList<MIX_Interval> intersec_al = new ArrayList<>();

                                            BufferedReader br_intersec = new BufferedReader(new FileReader(fileDirHM_Chr_sample));
                                            String line_intersec = "";
                                            while((line_intersec=br_intersec.readLine())!=null)
                                            {
                                                String[] split = line_intersec.split("\t");

                                                MIX_Interval_Object mio = new MIX_Interval_Object(split[0],Integer.parseInt(split[1]),Integer.parseInt(split[2]),split[3],Integer.parseInt(split[4]),split[5],Double.parseDouble(split[6]),Double.parseDouble(split[7]),Double.parseDouble(split[8]));
                                                MIX_Interval mi = new MIX_Interval(Integer.parseInt(split[1]),Integer.parseInt(split[2]));
                                                mi.merged_intervals.add(mio);

                                                intersec_al.add(mi);

                                            }
                                            br.close();


                                            intersec_vals.put(fileDirHM_Chr_sample.getName(),intersec_al);
                                        }*/

                                        }
                                    }

                                    //intersections sorting!
                                /*if(options_intern.mix_option.equals("INTERSECTION"))
                                {
                                    for(String k: intersec_vals.keySet())
                                    {
                                        Collections.sort(intersec_vals.get(k));
                                    }
                                }*/

                                    //unions sorting
                                    Collections.sort(all_intervals);

                                    Stack<MIX_Interval> stack_union = mergeIntervals(all_intervals);

                                    ArrayList<MIX_Interval> chr_unions = new ArrayList<>();

                                    int min_occurence = 0;
                                    if(options_intern.mix_occurence_intersection==-1)
                                    {
                                        min_occurence=sample_number;
                                    }
                                    else
                                    {
                                        if(sample_number<options_intern.mix_occurence_intersection)
                                        {
                                            min_occurence=sample_number;

                                        }
                                        else
                                        {
                                            min_occurence=options_intern.mix_occurence_intersection;
                                        }
                                    }

                                    while(!stack_union.isEmpty())
                                    {
                                        MIX_Interval t = stack_union.pop();
                                        t.calculate_mean("HM_LEVEL");

                                        if(options_intern.mix_option.equals("INTERSECTION"))
                                        {
                                            if(t.merged_intervals.size()>=min_occurence)
                                            {
                                                chr_unions.add(t);
                                            }
                                        }
                                        else
                                        {
                                            chr_unions.add(t);
                                        }
                                    }

                                    Collections.sort(chr_unions);

                                    all_chromosomes.put(chr,chr_unions);

                                    try {
                                        chr_alpha.add(Integer.parseInt(chr));
                                    }
                                    catch (Exception e)
                                    {
                                        chr_str.add(chr);
                                    }
                                }
                            }

                            Collections.sort(chr_alpha);
                            Collections.sort(chr_str);


                            //print Sample_unions => can be used if not HM_Level is used
                            BufferedWriter bw = new BufferedWriter(new FileWriter(f_output_union_samples_tp_hm.getAbsolutePath()+File.separator+timepoint+"_"+hm+"."+file_ending));

                            int peak_counter = 1;

                            for(int chr : chr_alpha)
                            {
                                ArrayList<MIX_Interval> x = all_chromosomes.get(""+chr);
                                for(int i = 0; i < x.size();i++)
                                {
                                    bw.write(x.get(i).meanToString(peak_counter));
                                    bw.newLine();
                                    peak_counter++;
                                }
                            }

                            for(String chr: chr_str)
                            {
                                ArrayList<MIX_Interval> x = all_chromosomes.get(chr);
                                for(int i = 0; i < x.size();i++)
                                {
                                    bw.write(x.get(i).meanToString(peak_counter));
                                    bw.newLine();
                                    peak_counter++;
                                }
                            }
                            bw.close();
                        }
                    }
                }
            }
        }
    }

    /**
     * read config file
     * @param check_options should options be checked for validity? Should be true if pipeline is run, should be false if analyse programms are run
     */
    public void read_config_file(boolean check_options) throws IOException {

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
                case "mix_level":
                    options_intern.mix_level=split[1].substring(1,split[1].length()-1);
                    break;
                case "mix_option":
                    options_intern.mix_option=split[1].substring(1,split[1].length()-1);
                    break;
                case "mix_occurence_intersection":
                    options_intern.mix_occurence_intersection=Integer.parseInt(split[1]);
                    break;
                case "black_list_dir":
                    options_intern.black_list_dir=split[1].substring(1,split[1].length()-1);
                    break;
                case "black_list_signals":
                    options_intern.black_list_signals=new HashSet<>(Arrays.asList(split[1].substring(1,split[1].length()-1).toUpperCase().split(";")));
                    break;
                case "deseq2_input_directory":
                    options_intern.deseq2_input_directory=split[1].substring(1,split[1].length()-1);
                    break;
                case "deseq2_input_gene_id":
                    options_intern.deseq2_input_gene_id=split[1].substring(1,split[1].length()-1);
                    break;
                case "deseq2_biomart_dataset_species":
                    options_intern.deseq2_biomart_dataset_species=split[1].substring(1,split[1].length()-1);
                    break;
                case "deseq2_count_threshold":
                    options_intern.deseq2_count_threshold=Integer.parseInt(split[1]);
                    break;
                case "tepic_input_directory":
                    options_intern.tepic_input_directory=split[1].substring(1,split[1].length()-1);
                    options_intern.tepic_input_original=split[1].substring(1,split[1].length()-1);
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
                case "tgen_consensus":
                    options_intern.tgen_consensus=Double.parseDouble(split[1]);
                    break;
                case "tgen_consensus_calc":
                    options_intern.tgen_consensus_calc=split[1].substring(1,split[1].length()-1);
                    break;
                case "tgen_no_closest_locus":
                    options_intern.tgen_no_closest_locus=Boolean.parseBoolean(split[1]);
                    break;
                case "tgen_self_regulatory":
                    options_intern.tgen_self_regulatory=Boolean.parseBoolean(split[1]);
                    break;
                case "tgen_no_closest_tss":
                    options_intern.tgen_no_closest_tss=Boolean.parseBoolean(split[1]);
                    break;
                case "tgen_max_link_distances":
                    options_intern.tgen_max_link_distances=Integer.parseInt(split[1]);
                    break;
                case "tgen_pvalue":
                    options_intern.tgen_pvalue=Double.parseDouble(split[1]);
                    break;
                case "tgen_mt_writing":
                    options_intern.tgen_mt_writing=split[1].substring(1,split[1].length()-1);
                    break;
                case "dynamite_preprocessing_integrate_data_geneIDs":
                    options_intern.dynamite_preprocessing_integrate_data_geneIDs=Integer.parseInt(split[1]);
                    break;
                case "dynamite_preprocessing_integrate_data_log2fc":
                    options_intern.dynamite_preprocessing_integrate_data_log2fc=Integer.parseInt(split[1]);
                    break;
                case "dynamite_preprocessing_integrate_data_consider_geneFile":
                    options_intern.dynamite_preprocessing_integrate_data_consider_geneFile=split[1].substring(1,split[1].length()-1);
                    break;
                case "dynamite_out_var":
                    options_intern.dynamite_out_var=split[1].substring(1,split[1].length()-1);
                    break;
                case "dynamite_cores":
                    options_intern.dynamite_cores=Integer.parseInt(split[1]);
                    break;
                case "dynamite_alpha":
                    options_intern.dynamite_alpha=Double.parseDouble(split[1]);
                    break;
                case "dynamite_testsize":
                    options_intern.dynamite_testsize=Double.parseDouble(split[1]);
                    break;
                case "dynamite_Ofolds":
                    options_intern.dynamite_Ofolds=Integer.parseInt(split[1]);
                    break;
                case "dynamite_Ifolds":
                    options_intern.dynamite_Ifolds=Integer.parseInt(split[1]);
                    break;
                case "dynamite_balanced":
                    options_intern.dynamite_balanced=Boolean.parseBoolean(split[1]);
                    break;
                case "dynamite_performance":
                    options_intern.dynamite_performance=Boolean.parseBoolean(split[1]);
                    break;
                case "dynamite_randomise":
                    options_intern.dynamite_randomise=Boolean.parseBoolean(split[1]);
                    break;
                case "plot_th_coefficient":
                    String[] split_coefficient_ths = split[1].split(";");
                    options_intern.plot_th_coefficient.clear();
                    for(String s: split_coefficient_ths)
                    {
                        options_intern.plot_th_coefficient.add(Double.parseDouble(s));
                    }
                    break;
                case "plot_cutoff_tps":
                    options_intern.plot_cutoff_tps=Integer.parseInt(split[1]);
                    break;
                case "plot_cutoff_hms":
                    options_intern.plot_cutoff_hms=Integer.parseInt(split[1]);
                    break;
                case "plot_cutoff_gcs":
                    options_intern.plot_cutoff_gcs=Integer.parseInt(split[1]);
                    break;
                case "plot_top_k_genes":
                    options_intern.plot_top_k_genes=Integer.parseInt(split[1]);
                    break;
                case "website_interesting_tfs":
                    String[] split_interesting_tfs = split[1].substring(1,split[1].length()-1).split(";");
                    options_intern.website_interesting_tfs.addAll(Arrays.asList(split_interesting_tfs));
                    break;
                case "html_report_interesting_tfs":
                    String[] split_interesting_tfs_2 = split[1].substring(1,split[1].length()-1).split(";");
                    options_intern.website_interesting_tfs.addAll(Arrays.asList(split_interesting_tfs_2));
                    break;
                default:
                    logger.logLine("Misformed cfg file - please use template of: /COM2POSE/config_templates/com2pose_template.cfg");
                    logger.logLine("Do not delete unused parameters in config data!");
                    System.exit(1);
            }
        }
        br.close();

        if(check_options)
        {
            boolean all_set = checkOptions();
            logger.logLine("Check config file parameters for validity");
            if(!all_set)
            {
                logger.logLine("Not all [REQ]uired options set. Please set them in config file");
                logger.logLine("Aborting COM2POSE");
                System.exit(1);
            }
            else
            {
                logger.logLine("Parameters in config file valid");
            }
        }


        logger.logLine("Reading config file finished - no errors detected");


    }

    private boolean checkOptions() throws IOException {

        boolean all_set = true;

        /**
         * mix histone options
         */

        if(!options_intern.mix_level.equals(""))
        {
            if(!options_intern.mix_level.equals("HM_LEVEL")&&!options_intern.mix_level.equals("SAMPLE_LEVEL"))
            {
                logger.logLine("[MIX]: mix_level parameter must be either HM_LEVEL or SAMPLE_LEVEL");
                all_set=false;
            }
            if(!options_intern.mix_option.equals("UNION")&&!options_intern.mix_option.equals("INTERSECTION"))
            {
                logger.logLine("[MIX]: mix_option parameter must be either UNION or INTERSECTION");
                all_set=false;
            }
        }

        /**
         * black list options
         */

        if(!options_intern.black_list_dir.equals(""))
        {
            File f = new File(options_intern.black_list_dir);
            if(!f.exists())
            {
                logger.logLine("[BLACKLIST] blacklist path does not exist!");
                all_set=false;
            }
            if(f.isDirectory())
            {
                logger.logLine("[BLACKLIST] blacklist path is a directory but must be a file");
                all_set=false;
            }

            if(options_intern.black_list_signals.isEmpty())
            {
                logger.logLine("[BLACKLIST] signals must not be empty");
                all_set=false;
            }
        }

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

        if(options_intern.deseq2_biomart_dataset_species.equals("") && options_intern.tepic_ensg_symbol.equals(""))
        {
            logger.logLine("[DESEQ2] deseq2_biomart_dataset_species must be filled if tepic_ensg_symbol is empty!");
            all_set=false;
        }

        if(!options_intern.deseq2_biomart_dataset_species.equals(""))
        {
            options_intern.tepic_ensg_symbol=options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_deseq2_preprocessing+File.separator+options_intern.file_suffix_deseq2_mapping;
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

        if(options_intern.tepic_gene_annot.equals(""))
        {
            logger.logLine("[TEPIC] Gene annotation path is not given");
            all_set=false;
        }
        else
        {
            File f = new File(options_intern.tepic_gene_annot);
            if(!f.exists())
            {
                logger.logLine("[TEPIC] Gene annotation path does not exist!");
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
        }

        //check for map ensg symbol
        if(options_intern.tepic_ensg_symbol.equals("") && options_intern.deseq2_biomart_dataset_species.equals(""))
        {
            logger.logLine("[TEPIC] No map of ENSG to Gene Symbol is given!");
        }
        else if(options_intern.deseq2_biomart_dataset_species.equals(""))
        {
            File f = new File(options_intern.tepic_ensg_symbol);
            if(!f.exists())
            {
                logger.logLine("[TEPIC] Map of ENSG to Gene Symbol file path does not exist!");
                all_set=false;
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
         * TGEN options
         */
        if(!options_intern.path_tgen.equals(""))
        {
            File file_tgen = new File(options_intern.path_tgen);
            if(!file_tgen.exists() || !file_tgen.isDirectory())
            {
                logger.logLine("[TGENE] TGene file directory does not exist or is not a directory!");
                all_set=false;
            }

            File tgene_dir = new File(options_intern.path_tgen+File.separator+"bin");

            if(!tgene_dir.exists())
            {
                logger.logLine("[TGENE] TGene binary directory cannot be found: " + tgene_dir.getAbsolutePath());
                all_set=false;
            }

            if(options_intern.tgen_mt_writing.equals(""))
            {
                logger.logLine("[TGENE] Please specify spelling of Mitochondrial DNA, e.g. M or MT (default: MT)");
                all_set=false;
            }

            if(options_intern.tgen_consensus==0.0)
            {
                logger.logLine("[TGENE] tgen_consensus must be in range ]0.0,1.0], it cannot be 0.0, if you do not want to use consensus set path_tgen=\"\"");
                all_set=false;
            }

            if(options_intern.tgen_self_regulatory)
            {
                if(options_intern.tgen_consensus_calc.equals(""))
                {
                    logger.logLine("[TGENE] tgen_consensus_calc must be set.");
                    all_set=false;
                }
            }
        }

        /**
         * DYNAMITE OPTIONS
         */

        if(!options_intern.dynamite_preprocessing_integrate_data_consider_geneFile.equals(""))
        {
            File f = new File(options_intern.dynamite_preprocessing_integrate_data_consider_geneFile);
            if(!f.exists())
            {
                logger.logLine("[DYNAMITE] preprocessing consider genes file for integrateData.py does not exist!");
                all_set=false;
            }
        }

        /**
         * PLOT OPTIONS
         */

        if(options_intern.plot_th_coefficient.isEmpty())
        {
            logger.logLine("[PLOTS] plot th coefficients is empty, please use at least one coefficient.");
            all_set=false;
        }
        if(options_intern.plot_cutoff_tps<1)
        {
            logger.logLine("[PLOTS] plot_cutoff_tps must be >= 1");
        }
        if(options_intern.plot_cutoff_hms<1)
        {
            logger.logLine("[PLOTS] plot_cutoff_hms must be >= 1");
        }
        if(options_intern.plot_cutoff_gcs<0)
        {
            logger.logLine("[PLOTS] plot_cutoff_gcs must be >= 0");
        }
        if(options_intern.plot_top_k_genes<1)
        {
            logger.logLine("[PLOTS] plot_top_k_genes must be >= 1");
        }

        return all_set;
    }

    private HashMap<String, HashMap<String, HashSet<String>>> checkGroupsTEPIC()
    {
        HashMap<String, HashMap<String, HashSet<String>>> groups = new HashMap<>();

        File folder = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tepic_output_raw);

        for(File fileDir : folder.listFiles())
        {
            if(fileDir.isDirectory())
            {
                HashMap<String,HashSet<String>> n_timepoint = new HashMap<>();
                String n_tp_name = fileDir.getName();
                for(File fileDir2 : fileDir.listFiles())
                {
                    String n_hm_name="";
                    HashSet<String> n_hm = new HashSet<>();
                    n_hm_name= fileDir2.getName();

                    for(File fileDir3 : fileDir2.listFiles())
                    {
                        //samples
                        if(fileDir3.isDirectory())
                        {
                            n_hm.add(fileDir3.getName());
                        }
                    }
                    n_timepoint.put(n_hm_name,n_hm);

                }
                groups.put(n_tp_name,n_timepoint);
            }
        }
        return groups;
    }

    /**
     * BLACKLIST preprocessing - creates binary tree friendly files of chromosomes
     */
    private ArrayList<BL_ranges_binary_tree> recursive_split_BL(ArrayList<BL_ranges_binary_tree> region, ArrayList<BL_ranges_binary_tree> newly_ordered)
    {
        if(region.size()>0)
        {
            BL_ranges_binary_tree median = region.get(region.size()/2);
            newly_ordered.add(median);

            ArrayList<BL_ranges_binary_tree> region_left=new ArrayList<>();

            for(int i = 0; i < region.size()/2; i++)
            {
                region_left.add(region.get(i));
            }
            ArrayList<BL_ranges_binary_tree> region_right=new ArrayList<>();
            for(int i = region.size()/2+1; i < region.size(); i++)
            {
                region_right.add(region.get(i));
            }

            recursive_split_BL(region_left,newly_ordered);
            recursive_split_BL(region_right,newly_ordered);

        }
        return newly_ordered;
    }

    /**
     * TGENE preprocessing - creates binary tree friendly files of chromosomes
     */
    private  ArrayList<ENSG_ranges_binary_trees> recursive_split(ArrayList<ENSG_ranges_binary_trees> region, ArrayList<ENSG_ranges_binary_trees> newly_ordered)
    {
        if(region.size()>0)
        {
            ENSG_ranges_binary_trees median = region.get(region.size()/2);
            newly_ordered.add(median);

            ArrayList<ENSG_ranges_binary_trees> region_left=new ArrayList<>();
            for(int i = 0; i < region.size()/2; i++)
            {
                region_left.add(region.get(i));
            }
            ArrayList<ENSG_ranges_binary_trees> region_right=new ArrayList<>();
            for(int i = region.size()/2+1; i < region.size(); i++)
            {
                region_right.add(region.get(i));
            }

            recursive_split(region_left,newly_ordered);
            recursive_split(region_right,newly_ordered);

        }
        return newly_ordered;
    }

    /**
     * mixes the samples of one folder into one file, based on mix_option (UNION or INTERSECTION)
     */
    private static Stack<MIX_Interval> mergeIntervals(ArrayList<MIX_Interval> interval)
    {
        Stack<MIX_Interval> stack = new Stack<>();

        if(stack.empty())
        {
            stack.push(interval.get(0));
        }

        for(int i = 1; i < interval.size(); i++)
        {
            MIX_Interval top = stack.peek();

            if(top.end < interval.get(i).start)
            {
                stack.push(interval.get(i));
            } else if(top.end < interval.get(i).end)
            {
                top.end = interval.get(i).end;
                top.merged_intervals.addAll(interval.get(i).merged_intervals);

                stack.pop();
                stack.push(top);
            }
            else
            {
                top.merged_intervals.addAll(interval.get(i).merged_intervals);

                stack.pop();
                stack.push(top);
            }
        }
        return stack;
    }

    private void write_target_genes_of_tf(File f_input_target_genes_hm_group_clash, String timepoint, File f_out_hm_th_file, String tf, HashMap<String,String>ensg_gene_symbol_map) throws IOException {

        File f_output = new File(f_out_hm_th_file+File.separator+timepoint);
        f_output.mkdir();

        if(!f_input_target_genes_hm_group_clash.exists())
        {
            String[] split = f_input_target_genes_hm_group_clash.getName().split("_");

            String[] split_slashes = f_input_target_genes_hm_group_clash.getAbsolutePath().split(File.separator);
            String path = "";
            for(int i = 0; i < split_slashes.length-1;i++)
            {
                path += File.separator+split_slashes[i];
            }
            path+=File.separator+split[1]+"_"+split[0];

            f_input_target_genes_hm_group_clash = new File(path);
        }

        File f_input = new File(f_input_target_genes_hm_group_clash.getAbsolutePath()+File.separator+timepoint);



        for(File fileDir: f_input.listFiles())
        {
            if(fileDir.isFile())
            {
                BufferedReader br = new BufferedReader(new FileReader(fileDir));
                String line = br.readLine();

                int interesting_column = -1;

                String[] header = line.split("\t");
                for(int i = 0; i < header.length;i++)
                {
                    if(header[i].toUpperCase().matches(".*"+tf.toUpperCase()+".*"))
                    {
                        interesting_column=i;
                    }
                }

                if(interesting_column==-1)
                {
                    return;
                }

                ArrayList<Gene_Affinity_Value> all_affinities = new ArrayList<>();

                while((line=br.readLine())!=null)
                {
                    String[] split = line.split("\t");

                    Gene_Affinity_Value gav = new Gene_Affinity_Value();
                    gav.gene_name = split[0];
                    if(ensg_gene_symbol_map.containsKey(gav.gene_name))
                    {
                        gav.gene_symbol=ensg_gene_symbol_map.get(gav.gene_name);
                    }
                    gav.affinity_value=Double.parseDouble(split[interesting_column]);
                    all_affinities.add(gav);
                }
                br.close();

                Collections.sort(all_affinities);

                BufferedWriter bw = new BufferedWriter(new FileWriter(new File(f_output.getAbsolutePath()+File.separator+tf+".csv")));

                bw.write("ENSG\tSYMBOL\tAFFINITY");
                bw.newLine();

                for(int i = 0; i < options_intern.plot_top_k_genes;i++)
                {
                    bw.write(all_affinities.get(i).toString());
                    bw.newLine();
                }

                bw.close();
            }
        }
    }

    private String write_table_html(Double d, String level) throws IOException {

        File input_dir_root = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_out_analysis_data+File.separator+options_intern.folder_out_analysis_data_WEBSITE_OVERVIEW);
        //HashMap<String,File> hm_file = new HashMap<>();
        HashMap<String,HashMap<String,Boolean>> hm_tf_found = new HashMap<>();

        String hm = "";

        for(File fileDir: input_dir_root.listFiles())
        {
            HashMap<String,Boolean> tf_found = new HashMap<>();

            File f = new File(fileDir.getAbsolutePath()+File.separator+d+File.separator+options_intern.file_suffix_website_analysis_tf_available);
            //hm_file.put(fileDir.getName(),f);
            BufferedReader br = new BufferedReader(new FileReader(f));
            String line = br.readLine();
            while((line=br.readLine())!=null)
            {
                String[] split = line.split("\t");
                tf_found.put(split[0],Boolean.parseBoolean(split[1]));
            }
            hm_tf_found.put(fileDir.getName(),tf_found);
            hm = fileDir.getName();


            br.close();
        }

        HashMap<String,Boolean> first = hm_tf_found.get(hm);

        StringBuilder sb = new StringBuilder();
        sb.append("<div class=\"w3-content\">");
        sb.append("<h2> Threshold:" + d+"</h2>\n");
        sb.append("\t\t<table style=\"width:80%\">\n");
        sb.append("\t\t\t<tr>\n");
        sb.append("\t\t\t\t<th>");
        sb.append("TF");
        sb.append("\t\t\t\t</th>\n");

        for(String hms_key : hm_tf_found.keySet())
        {
            sb.append("\t\t\t\t<th>");
            sb.append(hms_key);
            sb.append("\t\t\t\t</th>\n");

        }
        sb.append("\t\t\t</tr>\n");

        //File figures = new File(options_intern.path_to_COM2POSE+File.separator+"ext"+File.separator+"WEBSITE"+File.separator+"images");


        for(String tf : first.keySet())
        {
            sb.append("\t\t\t<tr>\n");

            sb.append("\t\t\t\t<th>");
            sb.append(tf);
            sb.append("\t\t\t\t</th>\n");

            for(String key_hm: hm_tf_found.keySet())
            {
                HashMap<String,Boolean> current_tf_list = hm_tf_found.get(key_hm);
                if(current_tf_list.get(tf))
                {
                    sb.append("\t\t\t\t<th>");
                    if(level.equals("HOME"))
                    {
                        sb.append("<img src=\""+options_intern.folder_out_website_basics+File.separator+options_intern.folder_out_website_basics_website+File.separator+options_intern.folder_out_website_basics_website_images+File.separator+"is_available.png"+"\" style=\"width:50px;height:50px;\"/>");

                    }
                    if(level.equals("THRESHOLD"))
                    {
                        sb.append("<img src=\".."+File.separator+".."+File.separator+options_intern.folder_out_website_basics+File.separator+options_intern.folder_out_website_basics_website+File.separator+options_intern.folder_out_website_basics_website_images+File.separator+"is_available.png"+"\" style=\"width:50px;height:50px;\"/>");

                    }
                    sb.append("\t\t\t\t</th>\n");
                }
                else
                {
                    sb.append("\t\t\t\t<th>");
                    if(level.equals("HOME"))
                    {
                        sb.append("<img src=\""+options_intern.folder_out_website_basics+File.separator+options_intern.folder_out_website_basics_website+File.separator+options_intern.folder_out_website_basics_website_images+File.separator+"not_available.png"+"\" style=\"width:50px;height:50px;\"/>");

                    }
                    if(level.equals("THRESHOLD"))
                    {
                        sb.append("<img src=\".."+File.separator+".."+File.separator+options_intern.folder_out_website_basics+File.separator+options_intern.folder_out_website_basics_website+File.separator+options_intern.folder_out_website_basics_website_images+File.separator+"not_available.png"+"\" style=\"width:50px;height:50px;\"/>");
                    }
                    sb.append("\t\t\t\t</th>\n");
                }
            }
            sb.append("\t\t\t</tr>\n");
        }

        sb.append("</table>");
        sb.append("<div>");

        return sb.toString();
    }

    private String get_header_html(String level, String which_analysis) {

        File f_website_css = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_out_website+File.separator+ options_intern.folder_out_website_basics+File.separator+options_intern.folder_out_website_basics_website+File.separator+options_intern.folder_out_website_basics_website_css);

        File f_output_website = new File(options_intern.com2pose_working_directory+ File.separator+options_intern.folder_out_website);
        File f_output_website_htmls = new File(f_output_website.getAbsolutePath()+File.separator+options_intern.folder_out_website_htmls_regression_coefficients);

        StringBuilder sb_home_front = new StringBuilder();
        sb_home_front.append("<!DOCTYPE html>\n" +
                "<html lang=\"en\">\n" +
                "<title>COM2POSE: "+which_analysis
                +"</title>\n" +
                "<meta charset=\"UTF-8\">\n" +
                "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n");
        sb_home_front.append("<head>\n<script>\n" +
                "  function expand_collapse(id, div_id) {\n" +
                "  var x = document.getElementById(id).getAttribute(\"aria-expanded\"); \n" +
                "  if (x == \"true\") \n" +
                "  {\n" +
                "  x = \"false\"\n" +
                "  } else {\n" +
                "  x = \"true\"\n" +
                "  }\n" +
                "  document.getElementById(id).setAttribute(\"aria-expanded\", x);\n" +

                "  var y = document.getElementById(div_id); \n" +
                "  if (y.style.display === \"none\") \n" +
                "  {\n" +
                "  y.style.display = \"block\"\n" +
                "  } else {\n" +
                "  y.style.display = \"none\"\n" +
                "  }\n" +
                "  }\n" +
                "</script>\n</head>\n");

        for(File fileDir:f_website_css.listFiles())
        {
            if(fileDir.isFile())
            {
                sb_home_front.append("<link rel=\"stylesheet\" href=\"");
                if(level.equals("HOME"))
                {
                    sb_home_front.append(options_intern.folder_out_website_basics+File.separator+options_intern.folder_out_website_basics_website+File.separator+options_intern.folder_out_website_basics_website_css+File.separator+fileDir.getName());
                }
                if(level.equals("THRESHOLD"))
                {
                    sb_home_front.append(".."+File.separator+".."+File.separator+options_intern.folder_out_website_basics+File.separator+options_intern.folder_out_website_basics_website+File.separator+options_intern.folder_out_website_basics_website_css+File.separator+fileDir.getName());

                }
                if(level.equals("TF"))
                {
                    sb_home_front.append(".."+File.separator+".."+File.separator+".."+File.separator+options_intern.folder_out_website_basics+File.separator+options_intern.folder_out_website_basics_website+File.separator+options_intern.folder_out_website_basics_website_css+File.separator+fileDir.getName());
                }
                sb_home_front.append("\">\n");
            }
        }

        String rel_path = "";
        if(level.equals("HOME"))
        {
            rel_path=options_intern.html_report_home_regression_coefficient_analysis;
        }
        if(level.equals("THRESHOLD"))
        {
            rel_path=".."+File.separator+".."+File.separator+options_intern.html_report_home_regression_coefficient_analysis;
        }
        if(level.equals("TF"))
        {
            rel_path=".."+File.separator+".."+File.separator+".."+File.separator+options_intern.html_report_home_regression_coefficient_analysis;
        }


        String rel_path_distribution_analysis = "";
        if(level.equals("HOME"))
        {
            rel_path_distribution_analysis=options_intern.html_report_home_regression_distribution_analysis;
        }
        if(level.equals("THRESHOLD"))
        {
            rel_path_distribution_analysis=".."+File.separator+".."+File.separator+options_intern.html_report_home_regression_distribution_analysis;
        }
        if(level.equals("TF"))
        {
            rel_path_distribution_analysis=".."+File.separator+".."+File.separator+".."+File.separator+options_intern.html_report_home_regression_distribution_analysis;
        }


        sb_home_front.append("<style>\n" +
                "body,h1,h2,h3,h4,h5,h6 {font-family: \"Lato\"; width: \"max-content\";}\n" +
                ".w3-bar,h1\n" +
                //,button {font-family: "Montserrat", sans-serif}
                ".fa-anchor,.fa-coffee {font-size:200px}\n" +
                ".button {\n" +
                "  background-color: #4CAF50; /* Green */\n" +
                "  border: none;\n" +
                "  color: white;\n" +
                "  padding: 14px 31px;\n" +
                "  text-align: center;\n" +
                "  text-decoration: none;\n" +
                "  display: inline-block;\n" +
                "  font-size: 15px;\n" +
                "  margin: 4px 2px;\n" +
                "  cursor: pointer;\n" +
                "  width: 130px;\n" +
                "}\n"+
                "</style>\n" +
                "<body>\n" +
                "  \n" +
                "<!-- Header -->\n" +
                "<header class=\"w3-container w3-red w3-center\" style=\"padding:128px 16px\">\n" +
                " <a href='"+rel_path+"' target='_blank'><button class=\"button\">HOME coefficient analysis</button></a>\n" +
                " <a href='"+rel_path_distribution_analysis+"' target='_blank'><button class=\"button\">HOME distribution analysis</button></a>\n");

        String rel_path_parameter = "";
        if(level.equals("HOME"))
        {
            rel_path_parameter="PARAMETERS.html";
        }
        if(level.equals("THRESHOLD"))
        {
            rel_path_parameter=".."+File.separator+".."+File.separator+"PARAMETERS.html";
        }
        if(level.equals("TF"))
        {
            rel_path_parameter=".."+File.separator+".."+File.separator+".."+File.separator+"PARAMETERS.html";
        }

        sb_home_front.append("<a href='"+rel_path_parameter+"' target='_blank'><button class=\"button\" id=\"button_parameters\" > Parameters</button></a>\n");

        sb_home_front.append(
                "  <h1 class=\"w3-margin w3-jumbo\">COM2POSE: <i>results overview</i></h1>\n" +
                "  \n" );

        if(which_analysis.equals(options_intern.analysis_types_distribution_analysis))
        {
            sb_home_front.append("  <h2 class=\"w3-center\" style='margin-left:27%'><i>TF-TG Score Distribution Analysis</i></h2>\n");
        }

        if(which_analysis.equals(options_intern.analysis_types_regression_coefficient_analysis))
        {
            sb_home_front.append("  <h2 class=\"w3-center\" style='margin-left:27%'><i>Regression Coefficient Analysis</i></h2>\n");

            sb_home_front.append("    <div class=\"dropdown\">\n" +
                    "  <h2>Please choose a threshold</h2>\n" +
                    "  <button class=\"dropbtn\">Threshold</button>\n" +
                    "  <div class=\"dropdown-content\">\n");



            for(Double d: options_intern.plot_th_coefficient)
            {
                File f_output_website_htmls_th = new File(f_output_website_htmls.getAbsolutePath()+File.separator+d);
                f_output_website_htmls_th.mkdir();

                if(!threshold_folders_filled)
                {
                    threshold_folders.add(f_output_website_htmls_th);
                }

                String rel_path_2 = "";
                if(level.equals("HOME"))
                {
                    rel_path_2=options_intern.folder_out_website_htmls_regression_coefficients+File.separator+d+File.separator+"threshold_"+d+"_overview.html";
                }
                if(level.equals("THRESHOLD"))
                {
                    rel_path_2=".."+File.separator+d+File.separator+"threshold_"+d+"_overview.html";

                }
                if(level.equals("TF"))
                {
                    rel_path_2=".."+File.separator+".."+File.separator+d+File.separator+"threshold_"+d+"_overview.html";
                }


                sb_home_front.append("<a href=\"");
                sb_home_front.append(rel_path_2+"\" target='_blank'>Coefficient " + d);
                sb_home_front.append("</a>\n");
            }
            sb_home_front.append("  </div>\n");

        }



        sb_home_front.append("</div>\n" +
                "</header>\n");

        return sb_home_front.toString();
    }

    private void write_python_script_distribution_analysis(File input_background_file, File input_tf_root, File output_plots, File output_script_file, File output_stats) throws IOException {
        HashMap<String,String> composed_tfs = new HashMap<>();
        HashMap<String,HashSet<String>> composed_tfs_tfs = new HashMap<>();

        BufferedReader br_composed_tfs = new BufferedReader(new FileReader(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tepic_postprocessing+File.separator+options_intern.folder_name_tepic_postprocessing_tfs+File.separator+options_intern.file_suffix_tepic_postprocessing_tfs_tfs));
        String line_composed_tfs = "";
        while((line_composed_tfs=br_composed_tfs.readLine())!=null)
        {
            String[] split = line_composed_tfs.split("\t");

            HashSet temp_set = new HashSet();
            for(int i = 1; i < split.length;i++)
            {
                composed_tfs.put(split[i],split[0]);
                temp_set.add(split[i]);
            }
            composed_tfs_tfs.put(split[0],temp_set);
        }
        br_composed_tfs.close();

        String imports = "import pip\n" +
                "\n" +
                "def import_or_install(package):\n" +
                "    try:\n" +
                "        __import__(package)\n" +
                "    except ImportError:\n" +
                "        pip.main(['install', package])\n" +
                "\n" +
                "import io\n" +
                "from base64 import b64encode\n" +
                "import_or_install(\"plotly.express\")\n" +
                "import plotly.express as px\n" +
                "import_or_install(\"dash\")\n" +
                "import_or_install(\"dash_core_components\")\n" +
                "import_or_install(\"dash_html_components\")\n" +
                "import dash_core_components as dcc\n" +
                "import dash_html_components as html\n" +
                "from dash.dependencies import Input, Output\n" +
                "import plotly.graph_objs as go\n" +
                "\n" +
                "import_or_install(\"pandas\")\n" +
                "import_or_install(\"seaborn\")\n" +
                "import_or_install(\"matplotlib.pyplot\")\n" +
                "import pandas as pd\n" +
                "import seaborn as sns\n" +
                "import matplotlib.pyplot as plt\n" +
                "sns.set_context(\"notebook\")\n" +
                "color = \"#A6CEE3\"\n" +
                "sns.set_context(\"talk\")\n" +
                "sns.set_style(\"whitegrid\")\n" +
                "\n" +
                "import_or_install(\"numpy\")\n" +
                "import numpy as np\n" +
                "import_or_install(\"sts\")\n" +
                "import statistics as sts\n"+
                "plt.figure(figsize=(20, 17))\n\n\n" +
                "df_interesting_stats=pd.DataFrame(columns=['label','sum_all_values','number_target_genes','mean','median','95_quantile','99_quantile'])\n" +
                "row_counter=0\n\n" +
                "";

        StringBuilder sb_all = new StringBuilder(imports);

        sb_all.append("background=pd.read_table('");
        sb_all.append(input_background_file);
        sb_all.append("', comment=\"#\", usecols=['TF_TG_SCORE']).sort_values(['TF_TG_SCORE'], ascending=False)\n");
        sb_all.append("background[\"label\"] = \"background\"\n" +
                "\n"+
                "background_sum = sum(background[\"TF_TG_SCORE\"])\n" +
                "background_length = len(background)\n" +
                "background_mean = background_sum/background_length\n" +
                "background_median=sts.median(background['TF_TG_SCORE'])\n"+
                "background_quantile=np.percentile(background[\"TF_TG_SCORE\"],95)\n"+
                "background_quantile_99=np.percentile(background[\"TF_TG_SCORE\"],99)\n" +
                "df_interesting_stats.loc[0]=['background',background_sum,background_length,background_mean,background_median,background_quantile,background_quantile_99]\n"+
                "row_counter=1"+
                "\n\n\n");

        for(File fileDir: input_tf_root.listFiles())
        {
            String name_tf = fileDir.getName().split("\\.")[0].split("_")[0];
            String name_composed ="";

            if(composed_tfs.containsKey(name_tf))
            {
                name_composed=composed_tfs.get(name_tf);
            }

            sb_all.append(name_tf);
            sb_all.append("=pd.read_table('");
            sb_all.append(fileDir.getAbsolutePath());
            sb_all.append("', comment=\"#\", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)\n");
            sb_all.append(name_tf);
            sb_all.append(".columns=['TF_TG_SCORE','label']\n");
            sb_all.append(name_tf+"_sum = sum("+name_tf+"['TF_TG_SCORE'])\n" +
                    name_tf+"_length = len("+name_tf+")\n" +
                    name_tf+"_mean = "+name_tf+"_sum/"+name_tf+"_length\n");
            sb_all.append(name_tf+"_quantile=np.percentile("+name_tf+"['TF_TG_SCORE'], 99)\n");
            sb_all.append(name_tf+"_quantile_95=np.percentile("+name_tf+"['TF_TG_SCORE'], 95)\n");
            sb_all.append(name_tf+"_median=sts.median("+name_tf+"['TF_TG_SCORE'])\n");
            sb_all.append("if("+name_tf+"_median > background_median and "+name_tf+"_quantile > background_quantile):\n");
            sb_all.append("    background_"+name_tf+" = pd.concat([background,"+name_tf+"],axis=0)\n" +
                    "    ax_"+name_tf+" = sns.boxplot(x=\"label\", y=\"TF_TG_SCORE\",data=background_"+name_tf+",palette=\"Set3\")\n" +
                    "    ax_"+name_tf+".set_yscale(\"log\")\n");
            if(name_composed.equals(""))
            {
                sb_all.append("    plt.savefig(f'"+output_plots.getAbsolutePath()+File.separator+name_tf+".png')\n");
            }
            else
            {
                sb_all.append("    plt.savefig(f'"+output_plots.getAbsolutePath()+File.separator+name_composed+".png')\n");
            }
            sb_all.append("    del background_"+name_tf+"\n"+
                    "    plt.clf()\n" +
                    "    row_counter=row_counter+1\n");
            if(name_composed.equals(""))
            {
                sb_all.append("    df_interesting_stats.loc[row_counter]=['"+name_tf+"',"+name_tf+"_sum,"+name_tf+"_length,"+name_tf+"_mean,"+name_tf+"_median,"+name_tf+"_quantile_95,"+name_tf+"_quantile]\n");
            }
            else
            {
                sb_all.append("    df_interesting_stats.loc[row_counter]=['"+name_composed+"',"+name_tf+"_sum,"+name_tf+"_length,"+name_tf+"_mean,"+name_tf+"_median,"+name_tf+"_quantile_95,"+name_tf+"_quantile]\n");
            }
            sb_all.append("del "+name_tf+"\n"
                    + "plt.figure(figsize=(20, 17))\n\n\n");
        }

        File f_stats_all_print = new File(output_stats.getAbsolutePath()+File.separator+options_intern.file_suffix_distribution_analysis_plot_stats);

        sb_all.append("df_interesting_stats.to_csv('"+f_stats_all_print.getAbsolutePath()+"',sep='\t')\n");

        BufferedWriter bw_all = new BufferedWriter(new FileWriter(output_script_file.getAbsolutePath()));
        bw_all.write(sb_all.toString());
        bw_all.close();
    }
}
