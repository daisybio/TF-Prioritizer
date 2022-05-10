package tfprio.wrapper;

import tfprio.ArgParser;
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

        ExecutionTimeMeasurement buildTimer = new ExecutionTimeMeasurement();
        executeAndWait("docker-compose -f " + f_compose.getAbsolutePath() + " build", logger, true);
        logger.info("Finished building container. Process took " + buildTimer.stopAndGetDeltaFormatted());

        Process container = execute("docker-compose -f " + f_compose.getAbsolutePath() + " up", logger, true);

        String containerID =
                new String(Runtime.getRuntime().exec("docker-compose -f " + f_compose.getAbsolutePath() + " ps -q").getInputStream()
                .readAllBytes()).replace("\n", "");
        File statsFile = extend(argParser.getWorkingDirectory(), "stats.tsv");
        File statsExecutable = extend(argParser.getWorkingDirectory(), ".logDockerStats.sh");
        String statLogCommand = readFile(configs.scriptTemplates.f_logDockerStats.get()).replace("{STATSFILE}",
                statsFile.getAbsolutePath()).replace("{CONTAINER_ID}", containerID);
        writeFile(statsExecutable, statLogCommand);
        Process statLogging = execute(statsExecutable, logger);

        logger.info("Container ID: " + containerID);
        logger.info("Stat logging PID: " + statLogging.pid());

        Runtime.getRuntime().addShutdownHook(new Thread(() -> {
            statLogging.destroy();
            container.destroy();
        }));

        container.waitFor();
    }

    private static File buildCompose(Configs configs, ArgParser argParser, Logger logger) throws IOException
    {
        StringBuilder sb_compose = new StringBuilder();

        String uid = new String(Runtime.getRuntime().exec("id -u " + System.getProperty("user.name")).getInputStream()
                .readAllBytes()).replace("\n", "");
        String gid = new String(Runtime.getRuntime().exec("id -g " + System.getProperty("user.name")).getInputStream()
                .readAllBytes()).replace("\n", "");

        sb_compose.append("version: \"2.2\"\n");
        sb_compose.append("services:\n");
        sb_compose.append("\tcom2pose:\n");
        sb_compose.append("\t\tbuild:\n");
        sb_compose.append("\t\t\tcontext: ").append(argParser.getSourceDirectory().getAbsolutePath()).append("\n");
        sb_compose.append("\t\t\targs:\n");
        sb_compose.append("\t\t\t\t- USER_ID=").append(uid).append("\n");
        sb_compose.append("\t\t\t\t- GROUP_ID=").append(gid).append("\n");

        int availableProcessors = Runtime.getRuntime().availableProcessors();
        com.sun.management.OperatingSystemMXBean os =
                (com.sun.management.OperatingSystemMXBean) java.lang.management.ManagementFactory.getOperatingSystemMXBean();
        int availableMemoryMb = (int) (os.getTotalPhysicalMemorySize() / 1e6);
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
        sb_compose.append("\t\t\t- _JAVA_OPTIONS=-Xmx").append(configs.general.memoryLimitMb.get()).append("m\n");

        sb_compose.append("\t\tvolumes:\n");

        for (InputFileStructure structure : configs.getInputFileStructure())
        {
            if (structure.isSet())
            {
                File f_docker = extend(configs.general.d_docker_input.get(), structure.get().getName());

                sb_compose.append("\t\t\t- ").append(structure.get().getAbsolutePath()).append(":")
                        .append(f_docker.getAbsolutePath()).append("\n");

                structure.setValue(f_docker);
            }
        }

        configs.tgene.pathToExecutable.setValue(new File("/srv/dependencies/meme"));
        configs.igv.pathToIGV.setValue(new File("/srv/dependencies/igv"));

        File f_compose = extend(argParser.getWorkingDirectory(), "docker-compose.yml");
        File f_configs = extend(argParser.getWorkingDirectory(), "dockerConfigs.json");
        configs.save(f_configs);

        sb_compose.append("\t\t\t- ").append(argParser.getWorkingDirectory().getAbsolutePath()).append(":")
                .append(configs.general.d_docker_wd.get().getAbsolutePath()).append("\n");

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
