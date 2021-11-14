package further_analysis_tools.TPM_GC_Filter_Analysis;

import com2pose.COM2POSE_lib;
import org.apache.commons.cli.*;
import util.Logger;
import util.Options_intern;

import java.io.*;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;

public class TPM_GC_Filter_Analysis {

    public static void main(String[] args) throws Exception {

        //need working_dir folder_names
        Options_intern options_intern = new Options_intern();

        parseArguments(args, options_intern);

        Logger logger = new Logger(options_intern.write_to_logfile, options_intern.tpm_gc_filter_analysis_working_dir);

        logger.logLine("#########################################");
        logger.logLine("############## COM2POSE #################");
        logger.logLine("#########################################");
        logger.logLine("######## TPM GC Filter analysis #########");
        logger.logLine("#########################################");
        logger.logLine("Working directory set to: " + options_intern.tpm_gc_filter_analysis_working_dir);
        logger.logLine("TF-list directory set to: " + options_intern.tpm_gc_filter_analysis_tf_list);

        logger.logLine("Check input parameters:");
        //check input
        File working_directory = new File(options_intern.tpm_gc_filter_analysis_working_dir);
        File tf_list = new File(options_intern.tpm_gc_filter_analysis_tf_list);

        check_parameters(working_directory, tf_list, logger);

        HashMap<String, COM2POSE_lib> file_com2pose_obj =
                find_prepare_config_files(logger, working_directory, options_intern);

        HashSet<String> interesting_tfs = get_interesting_tfs(logger, tf_list);

        File output_dir = new File(
                working_directory + File.separator + options_intern.tpm_gc_filter_analysis_folder_directory_name);
        output_dir.mkdir();
        File output_dir_data = new File(output_dir.getAbsolutePath() + File.separator +
                options_intern.tpm_gc_filter_analysis_folder_directory_name_data);
        output_dir_data.mkdir();
        File output_dir_rscripts = new File(output_dir.getAbsolutePath() + File.separator +
                options_intern.tpm_gc_filter_analysis_folder_directory_name_RScripts);
        output_dir_rscripts.mkdir();
        File output_dir_plots = new File(output_dir.getAbsolutePath() + File.separator +
                options_intern.tpm_gc_filter_analysis_folder_directory_name_plots);
        output_dir_plots.mkdir();

        create_tf_folders(logger, output_dir_data, interesting_tfs);
        create_tf_folders(logger, output_dir_rscripts, interesting_tfs);
        create_tf_folders(logger, output_dir_plots, interesting_tfs);


        create_data(logger, options_intern, interesting_tfs, file_com2pose_obj, output_dir_data);

        create_Rscripts_execute(logger, options_intern, output_dir_data, output_dir_rscripts, output_dir_plots);

        logger.logLine("Finished analyzing impact of filters on interesting TFs.");

    }

    private static void create_Rscripts_execute(Logger logger, Options_intern options_intern, File output_dir_data,
                                                File output_dir_rscripts, File output_dir_plots) throws Exception {
        logger.logLine("Start writing RScripts and execute them");

        //TFs
        for (File fileDir : output_dir_data.listFiles()) {
            if (fileDir.isDirectory()) {
                for (File fileDir_Group : fileDir.listFiles()) {
                    if (fileDir_Group.isDirectory()) {
                        File input = new File(fileDir_Group.getAbsolutePath() + File.separator +
                                options_intern.tpm_gc_filter_analysis_suffix_data);
                        if (!input.exists()) {
                            logger.logLine(
                                    "No data in " + fileDir.getName() + ": " + fileDir_Group.getName() + "! Skipping.");
                            continue;
                        }

                        File file_output = new File(
                                output_dir_rscripts.getAbsoluteFile() + File.separator + fileDir.getName() +
                                        File.separator + fileDir_Group.getName() + ".R");

                        StringBuilder sb = new StringBuilder();
                        sb.append("library(ggplot2)\n");
                        sb.append("data <- read.csv(\"");
                        sb.append(fileDir_Group.getAbsolutePath() + File.separator +
                                options_intern.tpm_gc_filter_analysis_suffix_data);
                        sb.append("\",sep = \"\\t\")\n");

                        sb.append("png(file=\"");
                        sb.append(output_dir_plots.getAbsolutePath() + File.separator + fileDir.getName() +
                                File.separator + fileDir_Group.getName() + ".png");
                        sb.append("\")\n");

                        sb.append("p <- ggplot(data, aes(factor(FILE), SCORE, fill = HM)) + \n");
                        sb.append("  geom_bar(stat=\"identity\", position = \"dodge\") + \n" +
                                "  scale_fill_brewer(palette = \"Set1\")\n");
                        sb.append("p  <- p + ggtitle(\"");
                        sb.append(fileDir.getName());
                        sb.append(" in ");
                        sb.append(fileDir_Group.getName());
                        sb.append("\")\n");
                        sb.append("p + theme(axis.text.x = element_text(angle = 90))");


                        BufferedWriter bw = new BufferedWriter(new FileWriter(file_output));
                        bw.write(sb.toString());
                        bw.close();

                        String command = "Rscript " + file_output;

                        logger.logLine("[Rscript] execute script: " + command);
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
        }
        logger.logLine("Finished writing RScripts and executing them!");
    }

    private static void create_data(Logger logger, Options_intern options_intern, HashSet<String> interesting_tfs,
                                    HashMap<String, COM2POSE_lib> file_com2pose_obj, File output_dir_data)
            throws IOException {
        logger.logLine("Create data for plotting");


        HashMap<String, HashMap<String, HashMap<String, HashMap<String, Double>>>> tf_group_hm_file_score =
                new HashMap<>();
        //create hashmap for all tfs
        for (String s : interesting_tfs) {
            HashMap<String, HashMap<String, HashMap<String, Double>>> x = new HashMap<>();
            tf_group_hm_file_score.put(s, x);
        }

        //put data into datastructure and create necessay folders
        for (String s : file_com2pose_obj.keySet()) {
            COM2POSE_lib com2pose = file_com2pose_obj.get(s);

            File raw_path_folder_input = new File(com2pose.options_intern.com2pose_working_directory + File.separator +
                    options_intern.folder_out_put_DYNAMITE);

            //HISTONE MODIFICATIONS
            for (File fileDir_HM : raw_path_folder_input.listFiles()) {
                if (fileDir_HM.isDirectory()) {
                    //Groups
                    for (File fileDir_HM_Group : fileDir_HM.listFiles()) {
                        if (fileDir_HM_Group.isDirectory()) {
                            File input_data_dir = new File(
                                    raw_path_folder_input.getAbsolutePath() + File.separator + fileDir_HM.getName() +
                                            File.separator + fileDir_HM_Group.getName());
                            File input_data = new File(input_data_dir.getAbsolutePath() + File.separator +
                                    options_intern.file_suffix_dynamite_output_to_be_plotted);

                            for (String ss : interesting_tfs) {
                                File f = new File(output_dir_data + File.separator + ss + File.separator +
                                        fileDir_HM_Group.getName());
                                f.mkdir();
                                HashMap<String, HashMap<String, HashMap<String, Double>>> x =
                                        tf_group_hm_file_score.get(ss);
                                HashMap<String, HashMap<String, Double>> x_group;

                                if (x.containsKey(fileDir_HM_Group.getName())) {
                                    x_group = x.get(fileDir_HM_Group.getName());

                                    HashMap<String, Double> x_hm;

                                    if (x_group.containsKey(com2pose.options_intern.com2pose_working_directory)) {
                                        x_hm = x_group.get(com2pose.options_intern.com2pose_working_directory);
                                    } else {
                                        x_hm = new HashMap<>();
                                        x_group.put(com2pose.options_intern.com2pose_working_directory, x_hm);
                                    }
                                } else {
                                    x_group = new HashMap<>();
                                    HashMap<String, Double> x_hm = new HashMap<>();
                                    x_group.put(com2pose.options_intern.com2pose_working_directory, x_hm);
                                    x.put(fileDir_HM_Group.getName(), x_group);
                                }
                            }


                            BufferedReader br = new BufferedReader(new FileReader(input_data));
                            double mean_sum_pos = 0.0;
                            double mean_sum_neg = 0.0;
                            int mean_count = 0;
                            String line = br.readLine();
                            while ((line = br.readLine()) != null) {
                                String[] split = line.split("\t");

                                String[] split_tf = split[0].split("_");

                                if (interesting_tfs.contains(split_tf[0].toUpperCase())) {
                                    try {
                                        tf_group_hm_file_score.get(split_tf[0].toUpperCase())
                                                .get(fileDir_HM_Group.getName())
                                                .get(com2pose.options_intern.com2pose_working_directory)
                                                .put(fileDir_HM.getName(), Double.parseDouble(split[1]));

                                        if (options_intern.tpm_gc_filter_analysis_count_zeros) {
                                            mean_count++;
                                        }
                                        double score = Double.parseDouble(split[1]);
                                        if (score < 0) {
                                            mean_sum_neg += score;
                                        } else {
                                            mean_sum_pos += score;
                                        }
                                    } catch (Exception e) {
                                        e.printStackTrace();
                                    }
                                }
                            }
                            if (mean_count == 0) {
                                mean_count++;
                            }
                            mean_sum_pos /= mean_count;
                            mean_sum_neg /= mean_count;

                            for (String ss : interesting_tfs) {
                                tf_group_hm_file_score.get(ss.toUpperCase()).get(fileDir_HM_Group.getName())
                                        .get(com2pose.options_intern.com2pose_working_directory)
                                        .put("MEAN_POS", mean_sum_pos);
                                tf_group_hm_file_score.get(ss.toUpperCase()).get(fileDir_HM_Group.getName())
                                        .get(com2pose.options_intern.com2pose_working_directory)
                                        .put("MEAN_NEG", mean_sum_neg);

                            }
                            br.close();
                        }
                    }
                }
            }
        }

        //write datastructure into R readable data
        for (String k_tf : tf_group_hm_file_score.keySet()) {
            HashMap<String, HashMap<String, HashMap<String, Double>>> group_file_hm = tf_group_hm_file_score.get(k_tf);

            for (String k_group : group_file_hm.keySet()) {
                HashMap<String, HashMap<String, Double>> file_hm = group_file_hm.get(k_group);

                StringBuilder sb = new StringBuilder();
                sb.append("FILE\tHM\tSCORE");
                sb.append("\n");

                for (String k_file : file_hm.keySet()) {
                    HashMap<String, Double> hm_score = file_hm.get(k_file);
                    String[] split_k_file = k_file.split(File.separator);
                    String file_name = split_k_file[split_k_file.length - 2];


                    for (String k_hm : hm_score.keySet()) {
                        sb.append(file_name);
                        sb.append("\t");
                        sb.append(k_hm);
                        sb.append("\t");
                        sb.append(hm_score.get(k_hm));
                        sb.append("\n");
                    }
                }

                File f_output = new File(
                        output_dir_data.getAbsolutePath() + File.separator + k_tf + File.separator + k_group +
                                File.separator + options_intern.tpm_gc_filter_analysis_suffix_data);
                BufferedWriter bw = new BufferedWriter(new FileWriter(f_output));
                bw.write(sb.toString());
                bw.close();
            }
        }

        logger.logLine("Finished creating data for plotting");
    }

    private static void create_tf_folders(Logger logger, File output_dir, HashSet<String> interesting_tfs)
            throws IOException {
        logger.logLine("Create TF folders for: " + output_dir.getName());

        for (String s : interesting_tfs) {
            File f = new File(output_dir.getAbsolutePath() + File.separator + s);
            f.mkdir();
        }

    }

    private static HashSet<String> get_interesting_tfs(Logger logger, File tf_list) throws IOException {
        HashSet<String> tfs = new HashSet<>();
        logger.logLine("Read TF-list file");

        BufferedReader br = new BufferedReader(new FileReader(tf_list));
        String line = br.readLine();
        while ((line = br.readLine()) != null) {
            tfs.add(line.toUpperCase());

        }
        logger.logLine("Finished reading TF-list file");
        return tfs;
    }

    private static HashMap<String, COM2POSE_lib> find_prepare_config_files(Logger logger, File working_directory,
                                                                           Options_intern options_intern)
            throws IOException {
        logger.logLine("Find config files for available parameters:");
        HashMap<String, COM2POSE_lib> file_com2pose_obj = new HashMap<>();
        for (File fileDir : working_directory.listFiles()) {

            if (fileDir.isDirectory() &&
                    !fileDir.getName().equals(options_intern.tpm_gc_filter_analysis_folder_directory_name)) {
                Options_intern options_this_run = new Options_intern();
                options_this_run.config_data_path =
                        fileDir.getAbsolutePath() + File.separator + "com2pose_template.cfg";
                options_this_run.com2pose_working_directory =
                        fileDir.getAbsolutePath() + File.separator + options_intern.folder_name_usual_working_dir_name;
                COM2POSE_lib com2pose = new COM2POSE_lib(options_this_run, logger);
                com2pose.read_config_file(false);

                file_com2pose_obj.put(fileDir.getAbsolutePath(), com2pose);
            }

        }

        for (String s : file_com2pose_obj.keySet()) {
            COM2POSE_lib com2pose_obj = file_com2pose_obj.get(s);
            logger.logLine("Found filter in " + s);
            logger.logLine("TPM: " + com2pose_obj.options_intern.tepic_tpm_cutoff + " and GC filter: " +
                    com2pose_obj.options_intern.deseq2_count_threshold);
        }

        return file_com2pose_obj;

    }

    private static void check_parameters(File working_directory, File tf_list, Logger logger) throws IOException {


        if (!working_directory.isDirectory()) {
            logger.logLine("root-run-directories must be a directory not a file!");
            System.exit(1);
        }
        if (!working_directory.exists()) {
            logger.logLine("root-run-directories does not exist!");
            System.exit(1);
        }
        if (tf_list.isDirectory()) {
            logger.logLine("TF-list must be a file not a directory!");
            System.exit(1);
        }
        if (!tf_list.exists()) {
            logger.logLine("TF-list does not exist!");
            System.exit(1);
        }
        logger.logLine("Input parameters ok!");
    }

    private static void parseArguments(String[] args, Options_intern options_intern) {
        Options options = new Options();

        Option opt_working_dir = new Option("r", "root-run-directories", true,
                "[REQ]: root-run-directories where all runs are in different folder including a working_dir folder in each folder");
        opt_working_dir.setRequired(true);
        options.addOption(opt_working_dir);

        Option opt_tf_list = new Option("t", "TF-list", true, "[REQ]: TF-list file directory");
        opt_tf_list.setRequired(true);
        options.addOption(opt_tf_list);

        Option opt_write_logfile = new Option("l", "write-log-file", false,
                "[OPT]: if flag is set no logfile will be written, default: logfile will be written");
        options.addOption(opt_write_logfile);

        Option opt_count_zeros = new Option("c", "count-zeros", false,
                "[OPT]: if flag is set TFs with a score of 0 will not be counted for means, default: TFs with zeros are counted");
        options.addOption(opt_count_zeros);

        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd;

        try {
            cmd = parser.parse(options, args);


            options_intern.tpm_gc_filter_analysis_working_dir = cmd.getOptionValue("root-run-directories");
            options_intern.tpm_gc_filter_analysis_tf_list = cmd.getOptionValue("TF-list");

            if (cmd.hasOption("write-log-file")) {
                options_intern.write_to_logfile = false;
            }

            if (cmd.hasOption("count-zeros")) {
                options_intern.tpm_gc_filter_analysis_count_zeros = false;
            }


        } catch (ParseException e) {
            System.out.println(e.getMessage());
            formatter.printHelp("-r <root-run-directories> -t <TF-list> [-l] [-c]", options);
            System.exit(1);
        }


    }


}