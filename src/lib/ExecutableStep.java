package lib;

import tfprio.TFPRIO;
import util.Configs.Config;
import util.ExecutionTimeMeasurement;
import util.Logger;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import static util.FileManagement.*;
import static util.Hashing.*;

public abstract class ExecutableStep
{
    protected ExecutorService executorService = Executors.newFixedThreadPool(TFPRIO.configs.general.threadLimit.get());

    private final int shutDownTimeOut = 5;
    private final TimeUnit shutDownTimeUnit = TimeUnit.MINUTES;
    protected final Logger logger = new Logger(this.getClass().getSimpleName());

    public boolean simulate()
    {
        logger.debug("Simulation starting.");
        updateInputDirectory();
        if (checkRequirements())
        {
            broadcastCreatedFileStructure();
            logger.debug("Simulation successful.");
            return true;
        } else
        {
            logger.error("Simulation failed.");
            return false;
        }
    }

    public void run()
    {
        ExecutionTimeMeasurement timer = new ExecutionTimeMeasurement();
        logger.info("Starting.");
        execute();
        shutdown();
        logger.info("Finished. Step took " + timer.stopAndGetDeltaSeconds() + " seconds.");
    }

    protected abstract Set<Config<File>> getRequiredFileStructure();

    protected abstract Set<Config<File>> getCreatedFileStructure();

    protected abstract Set<Config<?>> getRequiredConfigs();

    protected void updateInputDirectory()
    {
    }

    private void broadcastCreatedFileStructure()
    {
        TFPRIO.createdFileStructure.addAll(getCreatedFileStructure());
    }

    private boolean checkRequirements()
    {
        boolean allGood = true;
        Set<Config<?>> requiredConfigs = getRequiredConfigs();
        Set<Config<File>> requiredFileStructure = getRequiredFileStructure();

        for (Config<?> config : requiredConfigs)
        {
            if (!config.isSet())
            {
                allGood = false;
                logger.warn("A required config is not set: " + config.getName());
            }
        }

        for (Config<File> fileConfig : requiredFileStructure)
        {
            if (!TFPRIO.createdFileStructure.contains(fileConfig))
            {
                allGood = false;
                logger.warn("Required FileStructure has not been created: " + fileConfig.get().getAbsolutePath());
            }
        }
        return allGood;
    }

    private String hashInputs() throws IOException
    {
        Set<Config<File>> requiredFileStructure = getRequiredFileStructure();
        ArrayList<File> requiredFiles = new ArrayList<>();

        for (Config<File> fileConfig : requiredFileStructure)
        {
            requiredFiles.add(fileConfig.get());
        }

        return hashFiles(requiredFiles);
    }

    private String hashConfigs()
    {
        Set<Config<?>> requiredConfigs = getRequiredConfigs();
        return util.Hashing.hashConfigs(requiredConfigs);
    }

    private String hashOutputs() throws IOException
    {
        Set<Config<File>> createdFileStructure = getCreatedFileStructure();
        ArrayList<File> createdFiles = new ArrayList<>();

        for (Config<File> config : createdFileStructure)
        {
            createdFiles.add(config.get());
        }

        return hashFiles(createdFiles);
    }

    public boolean verifyHash()
    {
        try
        {
            File hashFile = getHashFile();
            if (!hashFile.exists() || !hashFile.canRead())
            {
                logger.warn("Cannot read hash file.");
                return false;
            }
            List<String> content = readLines(hashFile);
            assert content.size() == 3;

            String input = content.get(0);
            String config = content.get(1);
            String output = content.get(2);

            String inputHash = hashInputs();
            if (!inputHash.equals(input))
            {
                logger.warn("Input changed");
                return false;
            }

            String configHash = hashConfigs();
            if (!configHash.equals(config))
            {
                logger.warn("Configs changed");
                return false;
            }

            String outputHash = hashOutputs();
            if (!outputHash.equals(output))
            {
                logger.warn("Output has changed since last run");
                return false;
            }
            return true;
        } catch (IOException | AssertionError e)
        {
            e.printStackTrace();
        }
        return false;
    }

    public void createHash()
    {
        try
        {
            String inputHash = hashInputs();
            String configHash = hashConfigs();
            String outputHash = hashOutputs();

            File hashFile = getHashFile();
            makeSureFileExists(hashFile);
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(hashFile)))
            {
                writer.write(inputHash);
                writer.newLine();
                writer.write(configHash);
                writer.newLine();
                writer.write(outputHash);
            }
        } catch (IOException e)
        {
            e.printStackTrace();
        }
    }

    public void deleteHash()
    {
        File hashFile = getHashFile();
        if (hashFile.exists())
        {
            hashFile.delete();
        }
    }

    private File getHashFile()
    {
        return extend(TFPRIO.configs.general.d_workflowHashes.get(), this.getClass().getName() + ".md5");
    }

    protected void shutdown()
    {
        executorService.shutdown();

        try
        {
            executorService.awaitTermination(shutDownTimeOut, shutDownTimeUnit);
        } catch (InterruptedException e)
        {
            logger.error(
                    "Process did not finished within the defined limit of " + shutDownTimeOut + " " + shutDownTimeUnit);
            System.exit(1);
        }
    }

    protected void finishAllQueuedThreads()
    {
        shutdown();
        executorService = Executors.newFixedThreadPool(TFPRIO.configs.general.threadLimit.get());
    }

    protected abstract void execute();
}
