package util.Configs.Modules;

import org.json.JSONObject;
import tfprio.TFPRIO;
import util.Configs.Config;
import util.Logger;

import java.io.File;
import java.io.ObjectInputFilter;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.math.BigDecimal;
import java.util.HashMap;
import java.util.Map;

import static util.FileManagement.extend;


public abstract class AbstractModule
{
    protected final Config<File> workingDirectory;
    protected final Config<File> sourceDirectory;
    protected final Config<File> extDirectory;
    protected final Logger logger;

    protected Map<String, Config<?>> entries = new HashMap<>();
    protected Map<String, AbstractModule> subModules = new HashMap<>();

    public AbstractModule(Config<File> workingDirectory, Config<File> sourceDirectory, Logger logger)
    {
        this.workingDirectory = workingDirectory;
        this.sourceDirectory = sourceDirectory;
        this.extDirectory = extend(sourceDirectory, "ext");
        this.logger = logger;
    }

    protected void init()
            throws IllegalAccessException, InvocationTargetException, NoSuchMethodException, InstantiationException
    {
        registerEntries();
        initSubmodules();
    }

    private void initSubmodules()
            throws NoSuchMethodException, InvocationTargetException, InstantiationException, IllegalAccessException
    {
        Field[] fields = this.getClass().getFields();
        for (Field field : fields)
        {
            Class<?> superClass = field.getType().getSuperclass();

            if (superClass != null && superClass.equals(AbstractModule.class))
            {
                AbstractModule module =
                        (AbstractModule) field.getType().getConstructor(Config.class, Config.class, Logger.class)
                                .newInstance(workingDirectory, sourceDirectory, logger);
                field.set(this, module);
                subModules.put(field.getType().getSimpleName(), module);
            }
        }
    }

    private void registerEntries() throws IllegalAccessException
    {
        Field[] fields = this.getClass().getFields();
        for (Field field : fields)
        {
            if (field.getType().equals(Config.class))
            {
                entries.put(field.getName(), (Config<?>) field.get(this));
                ((Config<?>) field.get(this)).setName(this.getClass().getSimpleName() + " > " + field.getName());
            }
        }
    }

    public boolean merge(JSONObject mergeObject)
    {
        boolean worked = true;
        for (String key : mergeObject.keySet())
        {
            if (mergeObject.get(key).getClass().equals(JSONObject.class) && subModules.containsKey(key))
            {
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
                try
                {
                    entries.get(key).setValue(mergeObject.get(key));
                } catch (IllegalAccessException | ClassCastException | IllegalArgumentException e)
                {
                    worked = false;
                    logger.warn(this.getClass().getSimpleName() + ": " + key + ": " + e.getMessage());
                }
            } else
            {
                worked = false;
                logger.warn(this.getClass().getSimpleName() + ": Trying to set unknown config: " + key);
            }
        }
        return worked;
    }

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
                        !(fileConfig.get().getAbsolutePath().startsWith(workingDirectory.get().getAbsolutePath())))
                {
                    TFPRIO.createdFileStructure.add(fileConfig);
                }
            }
        }
        return subModulesValid && configsValid;
    }
}
