import com2pose.COM2POSE_lib;
import org.apache.commons.cli.*;


import util.Options_intern;

import java.io.File;
import java.io.IOException;

public class COM2POSE
{

    public static void main(String[] args) throws Exception {
        Options_intern options_intern= new Options_intern();
        parseArguments(args, options_intern);

        //prepare pipeline
        COM2POSE_lib com2pose_lib = new COM2POSE_lib(options_intern);
        com2pose_lib.read_config_file(true);

        //DESeq2
        //com2pose_lib.create_DESeq2_scripts();
        //com2pose_lib.run_and_postprocess_DESeq2();

        //TEPIC
        //com2pose_lib.run_tepic();
        //com2pose_lib.postprocess_tepic_output();

        //TGen
        if(!options_intern.path_tgen.equals(""))
        {
            //com2pose_lib.preprocess_tgen();
            //com2pose_lib.run_tgen();
            //com2pose_lib.merge_tgen();
            com2pose_lib.integrate_tgen();
        }


        //DYNAMITE
        //com2pose_lib.preprocess_dynamite();
        //com2pose_lib.run_DYNAMITE();

        //PLOTS
        //com2pose_lib.create_tp_plots();


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
