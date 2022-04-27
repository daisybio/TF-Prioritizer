package tfprio;

import org.apache.commons.cli.*;

import java.io.File;

public class ArgParser
{
    private File configFile;
    private File workingDirectory;
    private File sourceDirectory;

    public ArgParser(String[] args)
    {
        Options options = new Options();

        Option opt_configFile = new Option("c", "com2pose-config", true,
                "[REQ]: TFPRIO config file. Example in /TFPRIO/config_templates/configsTemplate.json");
        opt_configFile.setRequired(true);
        options.addOption(opt_configFile);

        Option opt_workingDirectory = new Option("w", "working-directory", true,
                "[REQ]: working directory where TFPRIO can create, remove and edit files");
        opt_workingDirectory.setRequired(true);
        options.addOption(opt_workingDirectory);

        Option opt_tfprioDirectory = new Option("p", "path-tfprio", true, "[REQ]: filepath to TFPRIO folder");
        opt_tfprioDirectory.setRequired(true);
        options.addOption(opt_tfprioDirectory);


        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd;

        try
        {
            cmd = parser.parse(options, args);

            configFile = new File(cmd.getOptionValue("com2pose-config"));
            workingDirectory = new File(cmd.getOptionValue("working-directory"));
            sourceDirectory = new File(cmd.getOptionValue("path-tfprio"));

        } catch (ParseException e)
        {
            System.out.println(e.getMessage());
            formatter.printHelp(
                    "-c <com2pose-config> -w <working-directory> -p <path-com2pose> [-t <tgen-dir>] [-l] [-m] [-a] [-t] [-b]",
                    options);
            System.exit(1);
        }
    }

    public File getConfigFile()
    {
        return configFile;
    }

    public File getWorkingDirectory()
    {
        return workingDirectory;
    }

    public File getSourceDirectory()
    {
        return sourceDirectory;
    }
}
