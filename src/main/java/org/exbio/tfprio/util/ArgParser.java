package org.exbio.tfprio.util;

import org.apache.commons.cli.*;

import java.io.File;

public class ArgParser {
    private final File configFile;
    private final File workingDirectory;

    public ArgParser(String[] args) throws ParseException {
        Options options = new Options();

        Option opt_configFile = new Option("c", "com2pose-config", true,
                "[REQ]: TFPRIO config file. Example in /TFPRIO/config_templates/configsTemplate.json");
        opt_configFile.setRequired(true);
        options.addOption(opt_configFile);

        Option opt_workingDirectory = new Option("w", "working-directory", true,
                "[REQ]: working directory where TFPRIO can create, remove and edit files");
        opt_workingDirectory.setRequired(true);
        options.addOption(opt_workingDirectory);

        CommandLineParser parser = new DefaultParser();
        CommandLine cmd;

        try {
            cmd = parser.parse(options, args);

            configFile = new File(cmd.getOptionValue("com2pose-config"));
            workingDirectory = new File(cmd.getOptionValue("working-directory"));
        } catch (ParseException e) {
            throw new ParseException("Failed to parse command line properties\n" + e.getMessage() + "\n" + options);
        }
    }

    public File getConfigFile() {
        return configFile;
    }

    public File getWorkingDirectory() {
        return workingDirectory;
    }
}
