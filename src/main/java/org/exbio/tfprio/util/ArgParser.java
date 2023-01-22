package org.exbio.tfprio.util;

import org.apache.commons.cli.*;

import java.io.File;

public class ArgParser {
    private final File configFile;
    private final File workingDirectory;
    private final Integer threadNumber;

    public ArgParser(String[] args) throws ParseException {
        Options options = new Options();

        Option opt_configFile = new Option("c", "config", true,
                "[REQ]: TFPRIO config file. Example in /TFPRIO/config_templates/configsTemplate.json");
        opt_configFile.setRequired(true);
        options.addOption(opt_configFile);

        Option opt_workingDirectory = new Option("o", "output-directory", true,
                "[REQ]: output directory where TF-PRIORITIZER can create, remove and edit files");
        opt_workingDirectory.setRequired(true);
        options.addOption(opt_workingDirectory);

        Option opt_threads = new Option("t", "thread-count", true, "[REQ]: thread count");
        opt_threads.setRequired(true);
        options.addOption(opt_threads);

        CommandLineParser parser = new DefaultParser();
        CommandLine cmd;

        try {
            cmd = parser.parse(options, args);

            configFile = new File(cmd.getOptionValue("config"));
            workingDirectory = new File(cmd.getOptionValue("output-directory"));
            threadNumber = Integer.parseInt(cmd.getOptionValue("thread-count"));
        } catch (ParseException e) {
            throw new ParseException("Failed to parse command line properties\n" + e.getMessage() + "\n" + options);
        }
    }

    public Integer getThreadNumber() {
        return threadNumber;
    }

    public File getConfigFile() {
        return configFile;
    }

    public File getWorkingDirectory() {
        return workingDirectory;
    }
}
