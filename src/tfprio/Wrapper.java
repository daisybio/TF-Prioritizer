package tfprio;

import util.Configs.ConfigTypes.InputFileStructure;
import util.Configs.Configs;
import util.Logger;

import java.io.File;
import java.io.IOException;

import static util.FileManagement.*;
import static util.ScriptExecution.executeAndWait;

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

        File f_compose = buildCompose(configs, argParser, logger);

        executeAndWait("docker-compose -f " + f_compose.getAbsolutePath() + " build", logger);
        executeAndWait("docker-compose -f " + f_compose.getAbsolutePath() + " up", logger);
    }

    private static File buildCompose(Configs configs, ArgParser argParser, Logger logger)
    {
        StringBuilder sb_compose = new StringBuilder();

        sb_compose.append("version: \"3\"\n");
        sb_compose.append("services:\n");
        sb_compose.append("\tcom2pose:\n");
        sb_compose.append("\t\tbuild: ").append(argParser.getSourceDirectory().getAbsolutePath()).append("\n");
        sb_compose.append("\t\tvolumes:\n");

        sb_compose.append("\t\t\t- ").append(argParser.getWorkingDirectory().getAbsolutePath()).append(":")
                .append(configs.general.d_docker_wd.get().getAbsolutePath()).append("\n");

        File d_input = extend(argParser.getWorkingDirectory(), "input");

        try
        {
            deleteFileStructure(d_input);
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        for (InputFileStructure structure : configs.getInputFileStructure())
        {
            if (structure.isSet())
            {
                File f_docker = extend(configs.general.d_docker_wd_input.get(), structure.get().getName());

                sb_compose.append("\t\t\t- ").append(structure.get().getAbsolutePath()).append(":")
                        .append(f_docker.getAbsolutePath()).append("\n");

                structure.setValue(f_docker);
            }
        }

        configs.tgene.pathToExecutable.setValue(new File("/srv/dependencies/meme"));
        configs.igv.pathToIGV.setValue(new File("/srv/dependencies/igv"));

        File f_compose = extend(argParser.getWorkingDirectory(), "docker-compose.yml");
        File f_configs = extend(d_input, "configs.json");
        configs.save(f_configs);

        try
        {
            writeFile(f_compose, sb_compose.toString().replace("\t", "  "));
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        return f_compose;
    }
}
