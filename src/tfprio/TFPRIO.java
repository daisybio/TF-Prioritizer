package tfprio;

import org.apache.commons.cli.*;

import util.Configs.Config;
import util.Configs.Configs;
import util.ExecutionTimeMeasurement;
import util.MapSymbolAndEnsg;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import static util.FileManagement.extend;


public class TFPRIO
{
    public static Configs configs;
    public static Map<String, Set<String>> groupsToHms = new HashMap<>();
    public static Map<String, Set<String>> groupCombinationsToHms = new HashMap<>();
    public static Set<Config<File>> createdFileStructure = new HashSet<>();
    public static MapSymbolAndEnsg mapSymbolAndEnsg;
    public static Config<File> latestInputDirectory;
    public static Set<String> existingHms = new HashSet<>();

    private static File workingDirectory;
    private static File tfprioDirectory;
    private static File configFile;

    public static void main(String[] args) throws Exception
    {
        ExecutionTimeMeasurement timer = new ExecutionTimeMeasurement();
        parseArguments(args);

        configs = new Configs(workingDirectory, tfprioDirectory);

        configs.merge(extend(tfprioDirectory, "config_templates", "configsTemplate.json"));
        configs.merge(configFile);
        configs.validate();

        Workflow workflow = new Workflow();

        if (workflow.simulationSuccessful())
        {
            workflow.run();
        }

        /*
        DESEQ2 Kann BatchVariablen akzeptieren
        Wenn entsprechende Config gesetzt ist, soll File übergeben werden, wo dann samples Batches zugeordnet sind.
        Diese Zuordnungen sollen dann DESeq2 übergeben werden
        */

        System.out.println("TFPRIO finished. Execution took " + timer.stopAndGetDeltaSeconds() + " seconds.");
    }


    private static void parseArguments(String[] args)
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

            TFPRIO.configFile = new File(cmd.getOptionValue("com2pose-config"));
            TFPRIO.workingDirectory = new File(cmd.getOptionValue("working-directory"));
            TFPRIO.tfprioDirectory = new File(cmd.getOptionValue("path-tfprio"));

        } catch (ParseException e)
        {
            System.out.println(e.getMessage());
            formatter.printHelp(
                    "-c <com2pose-config> -w <working-directory> -p <path-com2pose> [-t <tgen-dir>] [-l] [-m] [-a] [-t] [-b]",
                    options);
            System.exit(1);
        }
    }
}
