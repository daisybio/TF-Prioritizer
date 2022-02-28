package util.Configs.Modules;

import org.json.JSONObject;
import util.Configs.Config;
import util.Logger;

import java.io.File;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.util.HashMap;
import java.util.Map;

public abstract class AbstractModule
{
    protected final File workingDirectory;
    protected final File sourceDirectory;
    protected final Logger logger;

    protected Map<String, Config<?>> entries = new HashMap<>();
    protected Map<String, AbstractModule> subModules = new HashMap<>();

    public AbstractModule(File workingDirectory, File sourceDirectory, Logger logger)
    {
        this.workingDirectory = workingDirectory;
        this.sourceDirectory = sourceDirectory;
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
            Class superClass = field.getType().getSuperclass();

            if (superClass != null && superClass.equals(AbstractModule.class))
            {
                AbstractModule module =
                        (AbstractModule) field.getType().getConstructor(File.class, File.class, Logger.class)
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
            }
        }
    }

    public void merge(JSONObject mergeObject)
    {
        for (String key : mergeObject.keySet())
        {
            if (mergeObject.get(key).getClass().equals(JSONObject.class) && subModules.containsKey(key))
            {
                try
                {
                    subModules.get(key).merge(mergeObject.getJSONObject(key));
                } catch (ClassCastException e)
                {
                    logger.warn(key + ": " + e.getMessage());
                }

            } else if (entries.containsKey(key))
            {
                try
                {
                    entries.get(key).setValue(mergeObject.get(key));
                } catch (IllegalAccessException | ClassCastException e)
                {
                    logger.warn(key + ": " + e.getMessage());
                }
            } else
            {
                logger.warn("Trying to set unknown config: " + key);
            }
        }
    }

    public JSONObject toJSONObject()
    {

        JSONObject combined = new JSONObject();

        for (String key : subModules.keySet())
        {
            combined.accumulate(key, subModules.get(key).toJSONObject());
        }
        for (String key : entries.keySet())
        {
            combined.accumulate(key, entries.get(key).toString());
        }

        return combined;
    }

    protected Config<File> extend(Config<File> fileConfig, String extension)
    {
        return extend(fileConfig.get(), extension);
    }

    protected Config<File> extend(File file, String extension)
    {
        return new Config<>(new File(file.getAbsolutePath() + File.separator + extension));
    }
}
