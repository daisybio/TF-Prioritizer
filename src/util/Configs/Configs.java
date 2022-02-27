package util.Configs;

import util.Configs.Tools.AbstractTool;
import util.Configs.Tools.General;
import util.FileManagement;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.util.HashMap;
import java.util.Map;

import org.json.*;

public class Configs
{
    public final File workingDirectory, sourceDirectory;
    public final Map<String, AbstractTool> configs = new HashMap<>();
    public General general;

    public Configs(File workingDirectory, File sourceDirectory)
            throws ClassNotFoundException, NoSuchMethodException, IOException, InvocationTargetException,
            InstantiationException, IllegalAccessException
    {
        this.workingDirectory = workingDirectory;
        this.sourceDirectory = sourceDirectory;

        Field[] fields = this.getClass().getFields();
        for (Field field : fields)
        {
            Class superClass = field.getType().getSuperclass();

            if (superClass != null && superClass.equals(AbstractTool.class))
            {
                AbstractTool newTool = (AbstractTool) field.getType().getConstructor(File.class, File.class)
                        .newInstance(workingDirectory, sourceDirectory);
                field.set(this, newTool);
                configs.put(field.getType().getSimpleName(), newTool);
            }
        }
    }

    public void merge(File configFile) throws IOException
    {
        String content = FileManagement.loadFile(configFile);
        JSONObject allToolsJsonObject = new JSONObject(content);

        for (String toolName : allToolsJsonObject.keySet())
        {
            JSONObject toolJsonObject = allToolsJsonObject.getJSONObject(toolName);

            AbstractTool tool = configs.get(toolName);
            tool.merge(toolJsonObject);
        }
    }
}
