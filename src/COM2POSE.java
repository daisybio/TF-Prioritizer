import com2pose.COM2POSE_lib;
import org.apache.commons.cli.*;

import util.Options_intern;


public class COM2POSE
{

    public static void main(String[] args) throws Exception {
        Options_intern options_intern= new Options_intern();
        parseArguments(args, options_intern);


        //prepare pipeline
        COM2POSE_lib com2pose_lib = new COM2POSE_lib(options_intern);
        com2pose_lib.read_config_file(true);


        /*
        //mix histone modifications
        if(!options_intern.mix_level.equals(""))
        {
            com2pose_lib.mix_option();
        }


        //black listed regions filter
        if(!options_intern.black_list_dir.equals(""))
        {
            com2pose_lib.preprocess_blacklist();
            com2pose_lib.filter_blacklist();
        }

        if(options_intern.mix_mutually_exclusive)
        {
            com2pose_lib.mix_mutually_exclusive_peaks();
        }


        //DESeq2
        if(options_intern.tepic_ensg_symbol.equals("")||!options_intern.deseq2_biomart_dataset_species.equals(""))
        {
            com2pose_lib.get_ensg_symbol_mapping();
        }

        com2pose_lib.create_DESeq2_scripts();
        com2pose_lib.run_and_postprocess_DESeq2();


        if(!options_intern.path_tgen.equals("")) {
            com2pose_lib.preprocess_tgen();
            com2pose_lib.run_tgen();
            com2pose_lib.merge_tgen();
        }


        //TEPIC
        com2pose_lib.run_tepic();
        if(options_intern.tepic_randomize_tf_gene_matrix)
        {
            com2pose_lib.randomize_tepic();
        }
        com2pose_lib.postprocess_tepic_output();


        //TGen
        if(!options_intern.path_tgen.equals(""))
        {
            if(!options_intern.mix_mutually_exclusive)
            {
                com2pose_lib.create_tgen_groups();
            }
            com2pose_lib.filter_target_genes_tgen();

            if(options_intern.tgen_self_regulatory)
            {
                com2pose_lib.integrate_self_regulatory_tgen();
            }
        }

        //DYNAMITE
        com2pose_lib.preprocess_dynamite();
        com2pose_lib.run_DYNAMITE();

        //PLOTS
        com2pose_lib.create_tp_plots();
        com2pose_lib.analyze_plots_data();
        com2pose_lib.get_top_k_target_genes_plots();


        //CREATE OVERVIEW WEBSITE FOR EVERYTHING BEFORE DISTRIBUTION
        com2pose_lib.create_overview_html_report_before_distribution_analysis();

        //DISTRIBUTION ANALYSIS
        com2pose_lib.perform_distribution_analysis();
        com2pose_lib.create_distribution_plots();
        com2pose_lib.create_overview_html_report_distribution();*/
        //com2pose_lib.calculate_discounted_cumulative_gain_rank_distribution_analysis();

        if(!options_intern.igv_path_to_igv.equals(""))
        {
            com2pose_lib.run_igv();
        }

        System.out.println("X");
    }

    private static void parseArguments(String[] args, Options_intern options_intern)
    {
        Options options = new Options();

        Option opt_cfg_file= new Option("c","com2pose-config",true,"[REQ]: COM2POSE config file. Example in /COM2POSE/config_templates/com2pose_template.cfg");
        opt_cfg_file.setRequired(true);
        options.addOption(opt_cfg_file);

        Option opt_working_dir = new Option("w","working-directory",true,"[REQ]: working directory where COM2POSE can create, remove and edit files");
        opt_working_dir.setRequired(true);
        options.addOption(opt_working_dir);

        Option opt_path_to_COM2POSE = new Option("p","path-com2pose",true, "[REQ]: filepath to COM2POSE folder");
        opt_path_to_COM2POSE.setRequired(true);
        options.addOption(opt_path_to_COM2POSE);

        Option opt_tgen = new Option("t","tgen-dir", true, "[OPT]: use a consensus approach of TGen and TEPIC, provide directory to MEME suite installation folder here");
        options.addOption(opt_tgen);

        Option opt_write_logfile = new Option("l","write-log-file",false,"[OPT]: if flag is set no logfile will be written, default: logfile will be written");
        options.addOption(opt_write_logfile);


        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd;


        try
        {
            cmd = parser.parse(options, args);


            options_intern.config_data_path = cmd.getOptionValue("com2pose-config");
            options_intern.com2pose_working_directory=cmd.getOptionValue("working-directory");
            options_intern.path_to_COM2POSE=cmd.getOptionValue("path-com2pose");

            if(cmd.hasOption("tgen-dir"))
            {
                options_intern.path_tgen=cmd.getOptionValue("tgen-dir");
            }

            if(cmd.hasOption("write-log-file"))
            {
                options_intern.write_to_logfile=false;
            }

        } catch (ParseException e) {
            System.out.println(e.getMessage());
            formatter.printHelp("-c <com2pose-config> -w <working-directory> -p <path-com2pose> [-t <tgen-dir>] [-l]",options);
            System.exit(1);
        }
    }
}
