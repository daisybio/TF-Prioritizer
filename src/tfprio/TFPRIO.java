package tfprio;

import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.Configs;
import util.ExecutionTimeMeasurement;
import util.Logger;
import util.MapSymbolAndEnsg;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import static util.FileManagement.extend;


public class TFPRIO
{
    /**
     * The {@link Configs} object which is referenced by all the {@link lib.ExecutableStep} instances
     */
    public static Configs configs;

    /**
     * Mapping the existing groups to their associated histone modifications.
     */
    public static Map<String, Set<String>> groupsToHms = new HashMap<>();

    /**
     * Mapping the existing group combinations to their associated histone modifications.
     */
    public static Map<String, Set<String>> groupCombinationsToHms = new HashMap<>();

    public static Map<String, String> sample_group = new HashMap<>();
    public static Map<String, String> sample_batch = new HashMap<>();

    /**
     * Contains all the file configs whose allocated files or directories will be created during pipeline execution.
     * <p>
     * Is being filled during simulation.
     * In development mode, also existing file structures are considered as created by the pipeline. Otherwise,
     * commenting out executableSteps inside the {@link Workflow} would lead to a simulator error.
     */
    public static Set<AbstractConfig<File>> createdFileStructure = new HashSet<>();

    /**
     * Allows disk read efficient mapping of geneSymbols, geneIDs and their descriptions.
     */
    public static MapSymbolAndEnsg mapSymbolAndEnsg;

    /**
     * Contains all the histone modifications available in the input data.
     */
    public static Set<String> existingHms = new HashSet<>();

    public static File workingDirectory;
    public static File sourceDirectory;
    static File configFile;

    public static final boolean developmentMode = false;

    public static void main(String[] args) throws Exception
    {
        ExecutionTimeMeasurement timer = new ExecutionTimeMeasurement();
        ArgParser.parseArguments(args);

        configs = new Configs();

        configs.merge(extend(sourceDirectory, "config_templates", "defaultConfigs.json"));
        configs.merge(configFile);
        configs.validate();
        configs.save();

        Logger logger = new Logger("TFPRIO");

        logCurrentCommit(logger);

        if (developmentMode)
        {
            logger.warn("Development mode active.");
        }

        Workflow workflow = new Workflow();

        if (workflow.simulationSuccessful())
        {
            workflow.run();
        }

        logger.info("Finished. Execution took " + timer.stopAndGetDeltaFormatted() + ".");
    }

    private static void logCurrentCommit(Logger logger)
    {
        String command = "git rev-parse HEAD";

        try
        {
            Process p = Runtime.getRuntime().exec(command);
            BufferedReader input = new BufferedReader(new InputStreamReader(p.getInputStream()));

            String result = input.readLine();

            logger.info("Current commit: " + result);
        } catch (IOException e)
        {
            logger.warn("Could not get current commit: " + e.getMessage());
        }

    }
}
