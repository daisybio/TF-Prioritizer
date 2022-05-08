package tfprio;

import util.Configs.ConfigTypes.InputFileStructure;
import util.Configs.Configs;
import util.ExecutionTimeMeasurement;
import util.Logger;

import java.io.File;
import java.io.IOException;

import static util.FileManagement.*;
import static util.ScriptExecution.execute;
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

        File statsFile = extend(argParser.getWorkingDirectory(), "stats.tsv");
        File statsExecutable = extend(argParser.getWorkingDirectory(), ".logDockerStats.sh");
        String statLogCommand = readFile(configs.scriptTemplates.f_logDockerStats.get()).replace("{STATSFILE}",
                statsFile.getAbsolutePath());
        writeFile(statsExecutable, statLogCommand);

        ExecutionTimeMeasurement buildTimer = new ExecutionTimeMeasurement();
        executeAndWait("docker-compose -f " + f_compose.getAbsolutePath() + " build", logger, true);
        logger.info("Finished building container. Process took " + buildTimer.stopAndGetDeltaFormatted());
        Process statLogging = execute(statsExecutable, logger);
        executeAndWait("docker-compose -f " + f_compose.getAbsolutePath() + " up", logger, true);
        statLogging.destroy();
    }

    private static File buildCompose(Configs configs, ArgParser argParser, Logger logger)
    {
        StringBuilder sb_compose = new StringBuilder();

        sb_compose.append("version: \"2.2\"\n");
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

        int availableProcessors = Runtime.getRuntime().availableProcessors();
        com.sun.management.OperatingSystemMXBean os =
                (com.sun.management.OperatingSystemMXBean) java.lang.management.ManagementFactory.getOperatingSystemMXBean();
        int availableMemoryMb = (int) (os.getTotalMemorySize() / 1e6);
        int maxUsableMemoryMb = (int) (availableMemoryMb * 0.8);

        if (availableProcessors < configs.general.threadLimit.get())
        {
            logger.warn("The thread limit was set to " + configs.general.threadLimit.get() + " but this system only " +
                    "offers " + availableProcessors + " logical cores. The thread limit was reduced to " +
                    availableProcessors + ".");
            configs.general.threadLimit.setValue(availableProcessors);
        }
        if (maxUsableMemoryMb < configs.general.memoryLimitMb.get())
        {
            logger.warn("The memory limit was set to " + configs.general.memoryLimitMb.get() + "mb but this system " +
                    "only offers " + availableMemoryMb + "mb of memory. In order to leave some space for the os, the " +
                    "memory limit was set to " + maxUsableMemoryMb + "mb.");
            configs.general.memoryLimitMb.setValue(maxUsableMemoryMb);
        }

        if (configs.general.memoryLimitMb.get() < 8000)
        {
            logger.error("The memory limit has to be at least 8000mb. The set memory limit is " +
                    configs.general.memoryLimitMb + "mb.");
        }

        sb_compose.append("\t\tcpus: ").append(configs.general.threadLimit.get()).append("\n");
        sb_compose.append("\t\tmem_limit: ").append(configs.general.memoryLimitMb.get()).append("m\n");

        sb_compose.append("\t\tenvironment:\n");
        sb_compose.append("\t\t\t- _JAVA_OPTIONS=-Xmx").append(configs.general.memoryLimitMb.get()).append("m");

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
