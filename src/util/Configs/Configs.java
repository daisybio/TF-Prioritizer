package util.Configs;

import org.json.JSONException;
import org.json.JSONObject;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.Configs.ConfigTypes.InputFileStructure;
import util.Configs.ConfigTypes.SourceDirectoryFileStructure;
import util.Configs.Modules.AbstractModule;
import util.Configs.Modules.AngularReport.AngularReport;
import util.Configs.Modules.Blacklist.Blacklist;
import util.Configs.Modules.ChipAtlas.ChipAtlas;
import util.Configs.Modules.ConfusionMatrixes.ConfusionMatrixes;
import util.Configs.Modules.DeSeq2.DeSeq2;
import util.Configs.Modules.DistributionAnalysis.DistributionAnalysis;
import util.Configs.Modules.Dynamite.Dynamite;
import util.Configs.Modules.General;
import util.Configs.Modules.Igv.Igv;
import util.Configs.Modules.Jaspar.Jaspar;
import util.Configs.Modules.MixOptions.MixOptions;
import util.Configs.Modules.Plots.Plots;
import util.Configs.Modules.Report.Report;
import util.Configs.Modules.ScriptTemplates.ScriptTemplates;
import util.Configs.Modules.Tepic.Tepic;
import util.Configs.Modules.Tgene.Tgene;
import util.FileManagement;
import util.Logger;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import static util.FileManagement.*;

public class Configs
{
    /**
     * Maps the module names to their objects
     */
    private final Map<String, AbstractModule> configs = new HashMap<>();
    private final Logger logger;
    private final File logFile;
    /**
     * Configs that do not belong to a single module
     */
    public General general;
    /**
     * Configs for the Jaspar module
     */
    public Jaspar jaspar;
    /**
     * Configs for the Tgene module
     */
    public Tgene tgene;
    /**
     * Configs fot the Deseq2 module
     */
    public DeSeq2 deSeq2;
    /**
     * Configs for the Tepic module
     */
    public Tepic tepic;
    /**
     * Configs for the dynamite module
     */
    public Dynamite dynamite;
    /**
     * Configs for the IGV module
     */
    public Igv igv;
    /**
     * Configs for the Report module
     */
    public Report report;
    /**
     * Configs for the Blacklist module
     */
    public Blacklist blacklist;
    /**
     * Configs fot the MixOptions module
     */
    public MixOptions mixOptions;
    /**
     * Configs for the DistributionAnalysis module
     */
    public DistributionAnalysis distributionAnalysis;
    /**
     * Configs for the Plots module
     */
    public Plots plots;
    /**
     * Configs for the ScriptTemplates module
     */
    public ScriptTemplates scriptTemplates;
    /**
     * Configs for the ChipAtlas module
     */
    public ChipAtlas chipAtlas;
    /**
     * Configs for the AngularReport module
     */
    public AngularReport angularReport;
    public ConfusionMatrixes confusionMatrixes;

    /**
     * The default constructor.
     */
    public Configs(File workingDirectory, File sourceDirectory, File logFile)
    {
        File outputDirectory = extend(workingDirectory, "output");
        this.logFile = logFile;

        // Using the manual constructor since the TFPRIO.configs object is not yet available, since we are currently
        // inside the constructor of it.
        logger = new Logger("Configs", true, logFile);

        makeSureDirectoryExists(outputDirectory, logger);

        // Iterate all the fields inside this class
        Field[] fields = this.getClass().getFields();
        for (Field field : fields)
        {
            Class<?> superClass = field.getType().getSuperclass();

            // Check if the field extends AbstractModule
            if (superClass != null && superClass.equals(AbstractModule.class))
            {
                try
                {
                    // Call the AbstractModule constructor
                    AbstractModule module = (AbstractModule) field.getType()
                            .getConstructor(GeneratedFileStructure.class, SourceDirectoryFileStructure.class,
                                    Logger.class).newInstance(new GeneratedFileStructure(outputDirectory),
                                    new SourceDirectoryFileStructure(sourceDirectory), logger);
                    // Assign the created module to the field in this class
                    field.set(this, module);

                    // Store the created module in the name - object map
                    configs.put(field.getType().getSimpleName(), module);
                } catch (InvocationTargetException | InstantiationException | IllegalAccessException |
                         NoSuchMethodException e)
                {
                    e.printStackTrace();
                    // Problems during constructor call
                    logger.error(e.getMessage());
                }
            }
        }
    }

    /**
     * Merge an external config file to the configs object.
     *
     * @param configFile the config file
     * @throws IOException if the config file cannot be read
     */
    public void merge(File configFile) throws IOException
    {
        logger.debug("Merging configuration file: " + configFile.getAbsolutePath());
        String content = FileManagement.readFile(configFile);
        JSONObject combined = new JSONObject();
        boolean allModulesWorked = true;

        try
        {
            combined = new JSONObject(content);
        } catch (JSONException e)
        {
            logger.error("The config JSON-File does not match the JSON formant: " + e.getMessage());
        }

        for (String moduleName : combined.keySet())
        {
            JSONObject moduleJSONObject = combined.getJSONObject(moduleName);

            if (configs.containsKey(moduleName))
            {
                AbstractModule module = configs.get(moduleName);
                allModulesWorked = module.merge(moduleJSONObject) && allModulesWorked;
            } else
            {
                logger.warn("Trying to set config for unknown module: " + moduleName);
            }
        }
        if (!allModulesWorked)
        {
            logger.error("There were errors during config file merging. Aborting.");
        }
        logger.info("Merged configuration file: " + configFile.getAbsolutePath());
    }

    /**
     * Get a json string of all the stored configs (not only writable).
     *
     * @return the json string
     */
    public String toString()
    {
        return getConfigsJSONString(false);
    }

    /**
     * Get a json string of the stored configs.
     *
     * @param onlyWriteable defines if only writeable configs should be added to the config string
     * @return the json string
     */
    private String getConfigsJSONString(boolean onlyWriteable)
    {
        return getConfigsJSONObject(onlyWriteable).toString(4);
    }

    /**
     * Get a JSONObject of all the stored configs.
     *
     * @param onlyWriteable defines if only writeable configs should be added to the config string
     * @return the JSONObject
     */
    public JSONObject getConfigsJSONObject(boolean onlyWriteable)
    {
        return getConfigsJSONObject(onlyWriteable, false);
    }

    public JSONObject getConfigsJSONObject(boolean onlyWriteable, boolean excludeFiles)
    {
        JSONObject combined = new JSONObject();

        for (String key : configs.keySet())
        {
            JSONObject module = configs.get(key).toJSONObject(onlyWriteable, excludeFiles);
            if (!module.isEmpty())
            {
                combined.accumulate(key, module);
            }
        }

        return combined;
    }

    /**
     * Validate the configs inside all the modules.
     */
    public void validate()
    {
        logger.info("Validating configs");
        boolean allValid = true;
        for (AbstractModule module : configs.values())
        {
            allValid = module.validate() && allValid;
        }
        if (allValid)
        {
            logger.info("Configs are valid.");
        } else
        {
            logger.error("Configs are invalid. Aborting.");
        }
    }

    /**
     * Save the configs to a file in json format
     */
    public void save(File saveFile)
    {
        try
        {
            writeFile(saveFile, getConfigsJSONString(true));
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }
    }

    public File getLogFile()
    {
        return logFile;
    }

    public Set<InputFileStructure> getInputFileStructure()
    {
        Set<InputFileStructure> inputFileStructures = new HashSet<>();

        for (Map.Entry<String, AbstractModule> entry : configs.entrySet())
        {
            inputFileStructures.addAll(entry.getValue().getInputFileStructure());
        }

        return inputFileStructures;
    }
}
