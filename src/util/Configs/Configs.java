package util.Configs;

import util.Configs.Modules.*;
import util.Configs.Modules.Blacklist.Blacklist;
import util.Configs.Modules.ChipAtlas.ChipAtlas;
import util.Configs.Modules.DeSeq2.DeSeq2;
import util.Configs.Modules.DistributionAnalysis.DistributionAnalysis;
import util.Configs.Modules.Dynamite.Dynamite;
import util.Configs.Modules.Igv.Igv;
import util.Configs.Modules.Jaspar.Jaspar;
import util.Configs.Modules.MixOptions.MixOptions;
import util.Configs.Modules.Plots.Plots;
import util.Configs.Modules.Report.Report;
import util.Configs.Modules.ScriptTemplates.ScriptTemplates;
import util.Configs.Modules.Tepic.Tepic;
import util.Configs.Modules.Tgene.Tgene;
import util.Configs.Modules.TpmGcFilterAnalysis.TpmGcFilterAnalysis;
import util.FileManagement;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.util.HashMap;
import java.util.Map;

import org.json.*;
import util.Logger;

import static util.FileManagement.extend;

public class Configs
{
    public final Config<File> workingDirectory, sourceDirectory;
    private final Map<String, AbstractModule> configs = new HashMap<>();
    public General general;
    public Jaspar jaspar;
    public Tgene tgene;
    public DeSeq2 deSeq2;
    public Tepic tepic;
    public Dynamite dynamite;
    public Igv igv;
    public Report report;
    public Blacklist blacklist;
    public MixOptions mixOptions;
    public DistributionAnalysis distributionAnalysis;
    public Plots plots;
    public ScriptTemplates scriptTemplates;
    public ChipAtlas chipAtlas;
    public TpmGcFilterAnalysis tpmGcFilterAnalysis;
    public Misc misc;

    private final Logger logger;


    public Configs(File workingDirectory, File sourceDirectory)
            throws ClassNotFoundException, NoSuchMethodException, IOException, InvocationTargetException,
            InstantiationException, IllegalAccessException
    {
        this.workingDirectory = new Config<>(extend(workingDirectory, "working_dir"));
        this.sourceDirectory = new Config<>(sourceDirectory);

        logger = new Logger("Configs", true,
                new File(workingDirectory.getAbsolutePath() + File.separator + "logfile.txt"));

        Field[] fields = this.getClass().getFields();
        for (Field field : fields)
        {
            Class<?> superClass = field.getType().getSuperclass();

            if (superClass != null && superClass.equals(AbstractModule.class))
            {
                AbstractModule module =
                        (AbstractModule) field.getType().getConstructor(Config.class, Config.class, Logger.class)
                                .newInstance(new Config<>(workingDirectory), new Config<>(sourceDirectory), logger);
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

            if (configs.containsKey(moduleName))
            {
                AbstractModule module = configs.get(moduleName);
                module.merge(moduleJSONObject);
            } else
            {
                logger.warn("Trying to set config for unknown module: " + moduleName);
            }


        }
        logger.info("Merged configuration file from " + configFile.getAbsolutePath());
    }

    public String toString()
    {
        return getConfigsJSONString(false);
    }

    private String getConfigsJSONString(boolean onlyWriteable)
    {
        JSONObject combined = new JSONObject();

        for (String key : configs.keySet())
        {
            combined.accumulate(key, configs.get(key).toJSONObject(onlyWriteable));
        }
        return combined.toString(4);
    }

    public void save(File file, boolean onlyWriteable) throws IOException
    {
        FileManagement.writeFile(file, getConfigsJSONString(onlyWriteable));
        logger.info("Saved configuration JSON to " + file.getAbsolutePath());
    }
}
