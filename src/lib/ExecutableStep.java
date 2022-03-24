package lib;

import tfprio.TFPRIO;
import util.Configs.Config;
import util.ExecutionTimeMeasurement;
import util.Logger;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.HashSet;
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

    private final TimeUnit shutDownTimeUnit = TimeUnit.MINUTES;
    protected final Logger logger = new Logger(this.getClass().getSimpleName());

    public boolean simulate()
    {
        logger.debug("Simulation starting.");
        verifyConfigUsage();
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
        boolean executed = false;
        if (TFPRIO.configs.general.hashingEnabled.get())
        {
            logger.debug("Verifying hash...");
            if (!verifyHash())
            {
                logger.debug("Hash is invalid. Clearing output.");
                deleteAllOutputs();
                logger.debug("Cleared outputs. Starting full execution.");
                execute();
                executed = true;
            } else
            {
                logger.debug("Skipped execution since hash is valid.");
            }
        } else
        {
            execute();
            executed = true;
        }
        shutdown();

        if (TFPRIO.configs.general.hashingEnabled.get() && executed)
        {
            createHash();
        }
        logger.info("Finished. Step took " + timer.stopAndGetDeltaSeconds() + " seconds.");
    }

    protected abstract Set<Config<File>> getRequiredFileStructure();

    protected abstract Set<Config<File>> getCreatedFileStructure();

    protected abstract Set<Config<?>> getRequiredConfigs();

    protected Set<Config<?>> getOptionalConfigs()
    {
        return new HashSet<>();
    }

    protected void updateInputDirectory()
    {
    }

    private void broadcastCreatedFileStructure()
    {
        for (Config<File> createdStructure : getCreatedFileStructure())
        {
            if (TFPRIO.createdFileStructure.contains(createdStructure))
            {
                logger.error("Writing to already existing structure: " + createdStructure.getName());
            } else
            {
                TFPRIO.createdFileStructure.add(createdStructure);
            }
        }
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

    private void deleteAllOutputs()
    {
        Set<Config<File>> outputFiles = getCreatedFileStructure();
        for (Config<File> outputFile : outputFiles)
        {
            try
            {
                deleteFile(outputFile.get());
            } catch (IOException e)
            {
                logger.error(e.getMessage());
            }
        }
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

    private String hashRequiredConfigs()
    {
        Set<Config<?>> requiredConfigs = getRequiredConfigs();
        return util.Hashing.hashConfigs(requiredConfigs);
    }

    private String hashOptionalConfigs()
    {
        Set<Config<?>> optionalConfigs = getOptionalConfigs();
        return util.Hashing.hashConfigs(optionalConfigs);
    }

    public String hashOutputs() throws IOException
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
            assert content.size() == 4;

            String oldInputHash = content.get(0);
            String oldRequiredConfigHash = content.get(1);
            String oldOptionalConfigHash = content.get(2);
            String oldOutputHash = content.get(3);

            String inputHash = hashInputs();
            if (!inputHash.equals(oldInputHash))
            {
                logger.warn("Input changed");
                return false;
            }

            String configHash = hashRequiredConfigs();
            if (!configHash.equals(oldRequiredConfigHash))
            {
                logger.warn("Required configs changed");
                return false;
            }

            String optionalConfigHash = hashOptionalConfigs();
            if (!optionalConfigHash.equals(oldOptionalConfigHash))
            {
                logger.warn("Optional configs changed.");
                return false;
            }

            String outputHash = hashOutputs();
            if (!outputHash.equals(oldOutputHash))
            {
                logger.warn("Output has changed since last run.");
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
            String requiredConfigHash = hashRequiredConfigs();
            String optionalConfigHash = hashOptionalConfigs();
            String outputHash = hashOutputs();

            File hashFile = getHashFile();
            makeSureFileExists(hashFile);
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(hashFile)))
            {
                writer.write(inputHash);
                writer.newLine();
                writer.write(requiredConfigHash);
                writer.newLine();
                writer.write(optionalConfigHash);
                writer.newLine();
                writer.write(outputHash);
            }
        } catch (IOException e)
        {
            e.printStackTrace();
        }
    }

    private File getHashFile()
    {
        return extend(TFPRIO.configs.general.d_workflowHashes.get(), this.getClass().getName() + ".md5");
    }

    protected void shutdown()
    {
        executorService.shutdown();
        int timeOutMinutes = getShutDownTimeOutMinutes();

        try
        {
            boolean ret = executorService.awaitTermination(timeOutMinutes, shutDownTimeUnit);
            if (!ret)
            {
                logger.error("Process did not finished within the defined limit of " + timeOutMinutes + " " +
                        shutDownTimeUnit);
                System.exit(1);
            }
        } catch (InterruptedException e)
        {
            logger.error("Process was interrupted");
            System.exit(1);
        }
    }

    protected void finishAllQueuedThreads()
    {
        shutdown();
        executorService = Executors.newFixedThreadPool(TFPRIO.configs.general.threadLimit.get());
    }

    private void verifyConfigUsage()
    {
        Set<Config<?>> registeredConfigs = new HashSet<>();
        registeredConfigs.addAll(getCreatedFileStructure());
        registeredConfigs.addAll(getRequiredConfigs());
        registeredConfigs.addAll(getOptionalConfigs());
        registeredConfigs.addAll(getRequiredFileStructure());

        for (Field field : this.getClass().getDeclaredFields())
        {
            if (Config.class.isAssignableFrom(field.getType()))
            {
                try
                {
                    field.setAccessible(true);
                    Config<?> current = (Config<?>) field.get(this);
                    if (!registeredConfigs.contains(current))
                    {
                        logger.warn("Config not assigned to a category: " + field.getName());
                    }
                } catch (IllegalAccessException e)
                {
                    logger.error(e.getMessage());
                } finally
                {
                    field.setAccessible(false);
                }
            }
        }
    }

    protected int getShutDownTimeOutMinutes()
    {
        return 5;
    }

    protected abstract void execute();
}
