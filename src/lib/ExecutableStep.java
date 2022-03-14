package lib;

import com2pose.COM2POSE;
import util.Logger;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

public abstract class ExecutableStep
{
    protected ExecutorService executorService =
            Executors.newFixedThreadPool(COM2POSE.configs.general.threadLimit.get());

    private final int shutDownTimeOut = 5;
    private final TimeUnit shutDownTimeUnit = TimeUnit.MINUTES;
    protected final Logger logger = new Logger(this.getClass().getSimpleName());

    public ExecutableStep()
    {
        logger.info("Starting");
        long startTime = System.currentTimeMillis();
        execute();
        shutdown();
        double deltaSeconds = (double) (System.currentTimeMillis() - startTime) / 1e3;
        logger.info("Finished. Step took " + deltaSeconds + " seconds.");
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
        executorService = Executors.newFixedThreadPool(COM2POSE.configs.general.threadLimit.get());
    }

    protected abstract void execute();
}
