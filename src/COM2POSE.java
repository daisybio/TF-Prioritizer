import com2pose.COM2POSE_lib;
import org.apache.commons.cli.*;


import util.Options_intern;

import java.io.IOException;

public class COM2POSE
{

    public static void main(String[] args) throws IOException {
        Options_intern options_intern= new Options_intern();
        parseArguments(args, options_intern);

        COM2POSE_lib com2pose_lib = new COM2POSE_lib(options_intern);

        //DESeq2
        com2pose_lib.read_config_file();
        com2pose_lib.create_DESeq2_scripts();
        com2pose_lib.run_and_postprocess_DESeq2();

        //TEPIC

        //DYNAMITE

        //CREATE PLOTS





        System.out.println("X");
    }

    private static void parseArguments(String[] args, Options_intern options_intern)
    {
        Options options = new Options();

        Option opt_cfg_file= new Option("c","com2pose-config",true,"[REQ]: COM2POSE config file. Example in /COM2POSE/config_templates/com2pose_template.cfg");
        opt_cfg_file.setRequired(true);
        options.addOption(opt_cfg_file);

        /**
         * required options

        Option opt_input_dir_peaks = new Option("i", "input-dir-peaks", true, "[REQ] input path to directory with ChIP-seq peak data.");
        opt_input_dir_peaks.setRequired(true);
        options.addOption(opt_input_dir_peaks);

        Option opt_reference_genome = new Option("g", "input-reference-genome", true, "[REQ] filepath to reference genome");
        opt_reference_genome.setRequired(true);
        options.addOption(opt_reference_genome);

        Option opt_gene_annotation = new Option("a", "input-gene-annotation", true, "[REQ] filepath to gene annotation file");
        opt_gene_annotation.setRequired(true);
        options.addOption(opt_gene_annotation);

        Option opt_output_directory = new Option("o", "output-directory", true, "[REQ] directory for output files");
        opt_output_directory.setRequired(true);
        options.addOption(opt_output_directory);

        /**
         * optional options


        Option opt_window = new Option("w", "window-size", true, "[OPT] Size of the sliding window, default = 50000");
        options.addOption(opt_window);

        Option opt_cores = new Option("c","cores",true, "[OPT] how many cores should be used for parallelization, default = 1");
        options.addOption(opt_cores);

        /**
         * FLAGS

        Option opt_subdirs = new Option("s","compute-subdirectories",false,"[FLAG] compute all files in subdirectories, default = false");
        options.addOption(opt_subdirs);
        */

        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd;


        try
        {
            cmd = parser.parse(options, args);

            if(cmd.getOptionValue("com2pose-config")!=null)
            {
                options_intern.config_data_path = cmd.getOptionValue("com2pose-config");
            }
            else
            {
                System.err.println("COM2POSE config data must be set");
                System.exit(1);
            }

            /*
            //set intern options
            options_intern.input_dir_peaks = cmd.getOptionValue("input-dir-peaks");
            options_intern.input_reference_genome = cmd.getOptionValue("input-reference-genome");
            options_intern.output_dir = cmd.getOptionValue("output-directory");
            options_intern.gene_annotation = cmd.getOptionValue("input-gene-annotation");

            if(cmd.getOptionValue("window-size")!=null)
            {
                try {
                    options_intern.window = Integer.parseInt(cmd.getOptionValue("window-size"));
                }catch (Exception e)
                {
                    e.printStackTrace();
                    System.out.println("-w\t--window-size is not an Integer");
                    System.exit(1);
                }
            }

            if(cmd.getOptionValue("cores")!=null)
            {
                try{
                    options_intern.cores = Integer.parseInt(cmd.getOptionValue("cores"));
                }catch (Exception e)
                {
                    e.printStackTrace();
                    System.out.println("-c\t--cores is not an Integer");
                    System.exit(1);
                }
            }


            if(cmd.hasOption("compute-subdirectories"))
            {
                options_intern.run_all_subdirectories=true;

            }*/

        } catch (ParseException e) {
            System.out.println(e.getMessage());
            formatter.printHelp("-c <com2pose-config>",options);
            //formatter.printHelp("-i <input-dir-peaks> -g <input-reference-genome> -a <input-gene-annotation> -o <output-directory> [-w <window-size] [-c <cores>] [-s]", options);
            System.exit(1);
        }
    }
}
