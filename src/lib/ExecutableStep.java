package lib;

import tfprio.TFPRIO;
import util.Configs.Config;
import util.Logger;

import java.io.File;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

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
        logger.info("Starting.");
        long startTime = System.currentTimeMillis();
        execute();
        shutdown();
        double deltaSeconds = (double) (System.currentTimeMillis() - startTime) / 1e3;
        logger.info("Finished. Step took " + deltaSeconds + " seconds.");
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
