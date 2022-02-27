package util.Configs.Tools;

import org.json.JSONObject;
import util.Configs.Config;

import java.io.File;
import java.lang.reflect.Field;
import java.util.HashMap;
import java.util.Map;

public abstract class AbstractTool
{
    final File workingDirectory;
    final File sourceDirectory;

    Map<String, Config<?>> entries = new HashMap<>();

    public AbstractTool(File workingDirectory, File sourceDirectory)
    {
        this.workingDirectory = workingDirectory;
        this.sourceDirectory = sourceDirectory;
    }

    void registerEntries() throws IllegalAccessException
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
            if (entries.containsKey(key))
            {
                entries.get(key).setValue(mergeObject.get(key));
            }
        }
    }

    @Override public String toString()
    {
        return "AbstractTool{" + "entries=" + entries + '}';
    }
}
