package util.Configs;

import util.Configs.Modules.*;
import util.Configs.Modules.FileStructure.FileStructure;
import util.FileManagement;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.util.HashMap;
import java.util.Map;

import org.json.*;
import util.Logger;

public class Configs
{
    public final File workingDirectory, sourceDirectory;
    public final Map<String, AbstractModule> configs = new HashMap<>();
    public General general;
    public FileStructure fileStructure;
    public Jaspar jaspar;
    public Tgen tgen;
    public DeSeq2 deSeq2;
    public Tepic tepic;

    private final Logger logger;


    public Configs(File workingDirectory, File sourceDirectory)
            throws ClassNotFoundException, NoSuchMethodException, IOException, InvocationTargetException,
            InstantiationException, IllegalAccessException
    {
        this.workingDirectory = workingDirectory;
        this.sourceDirectory = sourceDirectory;
        logger = new Logger("Configs", true,
                new File(workingDirectory.getAbsolutePath() + File.separator + "logfile.txt"));

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
                configs.put(field.getType().getSimpleName(), module);
            }
        }
    }

    public void merge(File configFile) throws IOException
    {
        String content = FileManagement.loadFile(configFile);
        JSONObject combined = new JSONObject(content);

        for (String moduleName : combined.keySet())
        {
            JSONObject moduleJSONObject = combined.getJSONObject(moduleName);

            AbstractModule module = configs.get(moduleName);
            module.merge(moduleJSONObject);
        }
        logger.info("Merged configuration file from " + configFile.getAbsolutePath());
    }

    public String toString()
    {
        JSONObject combined = new JSONObject();

        for (String key : configs.keySet())
        {
            combined.accumulate(key, configs.get(key).toJSONObject());
        }
        return combined.toString(4);
    }

    public void save(File file) throws IOException
    {
        FileManagement.writeFile(file, toString());
        logger.info("Saved configuration JSON to " + file.getAbsolutePath());
    }
}
