package tfprio;

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

    static File workingDirectory;
    static File tfprioDirectory;
    static File configFile;

    public static void main(String[] args) throws Exception
    {
        ExecutionTimeMeasurement timer = new ExecutionTimeMeasurement();
        ArgParser.parseArguments(args);

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

        System.out.println("TFPRIO finished. Execution took " + timer.stopAndGetDeltaFormatted() + ".");
    }
}
