package util.Configs.Modules;

import org.json.JSONObject;
import tfprio.TFPRIO;
import util.Configs.Config;
import util.Logger;

import java.io.File;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.util.HashMap;
import java.util.Map;

import static util.FileManagement.extend;

/**
 * Abstract base class for config modules.
 * <p>
 * Contains configs and other modules (submodules).
 */
public abstract class AbstractModule
{
    protected final Config<File> workingDirectory;
    protected final Config<File> sourceDirectory;
    protected final Config<File> extDirectory;
    protected final Logger logger;

    /**
     * Maps the config names to their objects.
     * <p>
     * Required for json export.
     * The odd name is due to naming overlap with the {@link util.Configs.Configs} class.
     */
    protected Map<String, Config<?>> entries = new HashMap<>();

    /**
     * Maps all the submodule names inside this class to their objects.
     * <p>
     * Required for json export.
     */
    protected Map<String, AbstractModule> subModules = new HashMap<>();

    /**
     * The default constructor.
     *
     * @param workingDirectory the {@link TFPRIO} working directory
     * @param sourceDirectory  the {@link TFPRIO} source directory
     * @param logger           the {@link util.Configs.Configs} logger
     */
    public AbstractModule(Config<File> workingDirectory, Config<File> sourceDirectory, Logger logger)
    {
        this.workingDirectory = workingDirectory;
        this.sourceDirectory = sourceDirectory;
        this.extDirectory = extend(sourceDirectory, "ext");
        this.logger = logger;
    }

    /**
     * Fills the entry and submodule maps and does the same for all submodules.
     */
    protected void init()
            throws IllegalAccessException, InvocationTargetException, NoSuchMethodException, InstantiationException
    {
        registerEntries();
        initSubmodules();
    }

    /**
     * Calls the default constructor of all the submodules and registers them inside the submodule map.
     */
    private void initSubmodules()
            throws NoSuchMethodException, InvocationTargetException, InstantiationException, IllegalAccessException
    {
        // Get all fields inside the class extending the AbstractModule class
        Field[] fields = this.getClass().getFields();
        for (Field field : fields)
        {
            Class<?> superClass = field.getType().getSuperclass();

            // If the field superclass is AbstractModule
            if (superClass != null && superClass.equals(AbstractModule.class))
            {
                // Call the default constructor with the same argument as this object has been created
                AbstractModule module =
                        (AbstractModule) field.getType().getConstructor(Config.class, Config.class, Logger.class)
                                .newInstance(workingDirectory, sourceDirectory, logger);
                // Assign the created object to this object
                field.set(this, module);

                // Add the created object to the submodule map
                subModules.put(field.getType().getSimpleName(), module);
            }
        }
    }

    /**
     * Add all created configs to the entry map.
     */
    private void registerEntries() throws IllegalAccessException
    {
        // Get all fields of the class extending AbstractModule
        Field[] fields = this.getClass().getFields();
        for (Field field : fields)
        {
            // Check if the field is a Config
            if (field.getType().equals(Config.class))
            {
                // Add the config to the entry map
                entries.put(field.getName(), (Config<?>) field.get(this));

                // Store the config name inside the config. Makes logging easier.
                ((Config<?>) field.get(this)).setName(this.getClass().getSimpleName() + " > " + field.getName());
            }
        }
    }

    /**
     * Merge a config json object into the module.
     *
     * @param mergeObject the json object to merge
     * @return true if merging is successful, otherwise false
     */
    public boolean merge(JSONObject mergeObject)
    {
        boolean worked = true;
        for (String key : mergeObject.keySet())
        {
            if (mergeObject.get(key).getClass().equals(JSONObject.class) && subModules.containsKey(key))
            {
                // Merge configs to a submodule
                try
                {
                    subModules.get(key).merge(mergeObject.getJSONObject(key));
                } catch (ClassCastException e)
                {
                    worked = false;
                    logger.warn(this.getClass().getSimpleName() + ": " + key + ": " + e.getMessage());
                }

            } else if (entries.containsKey(key))
            {
                // Set config value
                try
                {
                    entries.get(key).setValueObject(mergeObject.get(key));
                } catch (IllegalAccessException | ClassCastException | IllegalArgumentException e)
                {
                    worked = false;
                    logger.warn(this.getClass().getSimpleName() + ": " + key + ": " + e.getMessage());
                }
            } else
            {
                // Key is neither a submodule nor an entry
                worked = false;
                logger.warn(this.getClass().getSimpleName() + ": Trying to set unknown config: " + key);
            }
        }
        return worked;
    }

    /**
     * Creates a json object containing the configs inside this instance and all submodules.
     *
     * @param onlyWriteable defines if only writeable configs should be included
     * @return the json object
     */
    public JSONObject toJSONObject(boolean onlyWriteable)
    {
        JSONObject combined = new JSONObject();

        for (String key : subModules.keySet())
        {
            JSONObject subModuleJSONObject = subModules.get(key).toJSONObject(onlyWriteable);
            if (!subModuleJSONObject.keySet().isEmpty())
            {
                combined.accumulate(key, subModuleJSONObject);
            }
        }
        for (String key : entries.keySet())
        {
            Config<?> entry = entries.get(key);

            if (onlyWriteable && !entry.isWriteable())
            {
                continue;
            }

            combined.accumulate(key, entry.toJSONifyAble());
        }

        return combined;
    }

    /**
     * Checks if the configs inside this module are valid. Does the same for all submodules.
     * <p>
     * If development mode is active and file configs are sub file structures of the working directory, they are
     * added to the {@link TFPRIO} createdFileStructure. This allows commenting out executableSteps from the
     * workspace without running into errors, if the files have already been created in an earlier pipeline execution
     * or are copied from another run
     *
     * @return true if all configs are valid, otherwise false
     */
    public boolean validate()
    {
        boolean subModulesValid = true;
        for (AbstractModule subModule : subModules.values())
        {
            subModulesValid = subModule.validate() && subModulesValid;
        }
        boolean configsValid = true;
        for (Config<?> config : entries.values())
        {
            configsValid = config.isValid() && configsValid;

            if (config.isValid() && config.isSet() && config.get().getClass().equals(File.class))
            {
                Config<File> fileConfig = (Config<File>) config;

                if (fileConfig.get().exists() &&
                        (!(fileConfig.get().getAbsolutePath().startsWith(workingDirectory.get().getAbsolutePath())) ||
                                TFPRIO.configs.general.developmentMode.get()))
                {
                    TFPRIO.createdFileStructure.add(fileConfig);
                }
            }
        }
        return subModulesValid && configsValid;
    }
}
