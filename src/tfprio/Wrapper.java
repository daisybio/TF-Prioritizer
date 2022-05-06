package tfprio;

import util.Configs.ConfigTypes.InputFileStructure;
import util.Configs.Configs;
import util.Logger;

import java.io.File;
import java.io.IOException;

import static util.FileManagement.*;

public class Wrapper
{
    public static void main(String[] args) throws Exception
    {
        ArgParser argParser = new ArgParser(args);

        File logFile = extend(argParser.getWorkingDirectory(), "wrapLog.txt");

        Logger logger = new Logger("Wrapper", true, logFile);

        Configs configs = new Configs(argParser.getWorkingDirectory(), argParser.getSourceDirectory(), logFile);
        try
        {
            configs.merge(extend(argParser.getSourceDirectory(), "config_templates", "defaultConfigs.json"));
            configs.merge(argParser.getConfigFile());
        } catch (IOException e)
        {
            logger.error("Exception during config merging: " + e.getMessage());
        }
        configs.validate();

        File d_links = extend(argParser.getWorkingDirectory(), "input");

        try
        {
            deleteFileStructure(d_links);
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        for (InputFileStructure structure : configs.getInputFileStructure())
        {
            if (structure.isSet())
            {
                File f_link = extend(d_links, structure.get().getName());
                softLink(f_link, structure.get(), logger);
                structure.setValue(f_link);
            }
        }
        configs.save(extend(d_links, "configs.json"));

        TFPRIO.execute(argParser);
    }
}
