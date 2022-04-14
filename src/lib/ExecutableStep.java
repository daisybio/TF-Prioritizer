package lib;

import tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
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
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;

import static util.FileManagement.*;
import static util.Hashing.*;

/**
 * Abstract class for a single analysis performed as part of the pipeline.
 * <p>
 * The following practices should be applied to all extending classes:
 * <ul>
 *     <li>All the used configs should be stored as private final class data elements in the beginning of the class</li>
 *     <li>The configs should be split to four blocks:
 *     <ol>
 *         <li>Required file structure (input files/directories)</li>
 *         <li>Created file structure (output files/directories)</li>
 *         <li>Required configs (mandatory configs for this executableStep)</li>
 *         <li>Optional configs (not mandatory but influencing the output)</li>
 *     </ol>
 *     </li>
 *     <li>Each config has to be assigned to one of the following methods:
 *     <ul>
 *         <li>{@link #getRequiredFileStructure}</li>
 *         <li>{@link #getCreatedFileStructure}</li>
 *         <li>{@link #getRequiredConfigs}</li>
 *         <li>{@link #getOptionalConfigs}</li>
 *     </ul>
 *     </li>
 *     <li>If a file structure or config is not required or created if a certain config constellation is active,
 *     the config constellation should be modelled inside the corresponding get method. The get methods should model
 *     the real execution requirements and outputs as exact as possible.</li>
 * </ul>
 */
public abstract class ExecutableStep
{
    /**
     * Allows multithreading with a defined number of threads.
     * The thread number can be set by overriding the {@link #getThreadNumber} method.
     */
    protected ThreadPoolExecutor executorService = (ThreadPoolExecutor) Executors.newFixedThreadPool(getThreadNumber());

    /**
     * The logger of this ExecutableStep.
     */
    protected final Logger logger = new Logger(this.getClass().getName().replace("lib.", "").replace("tfprio.", ""));

    private String noExecutionReason = null;

    /**
     * Check if all the requirements of this executableStep are met and broadcast the created file structures. Does
     * not execute the executableStep.
     *
     * @return true if the simulation was successful, otherwise false.
     */
    public boolean simulate()
    {
        logger.debug("Simulation starting.");
        updateInputDirectory();

        if (TFPRIO.developmentMode)
        {
            verifyConfigUsage();
        }

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

    /**
     * Wraps the executableStep execution with some framework checks.
     * <p>
     * Skips the executableStep if developmentMode is not active and valid hashes are found.
     * Stores new hashes if the executableStep has been executed and developmentMode is disabled.
     */
    public void run()
    {
        ExecutionTimeMeasurement timer = new ExecutionTimeMeasurement();
        logger.info("Starting.");
        boolean executed;
        if (!TFPRIO.developmentMode)
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
                executed = false;
            }
        } else
        {
            execute();
            executed = true;
        }
        shutdown();

        if (executed && !TFPRIO.developmentMode)
        {
            logger.debug("Writing hash");
            createHash();
        }
        logger.info("Finished. Step took " + timer.stopAndGetDeltaFormatted());
    }

    /**
     * Get the file configs whose value files or directories must have been created before execution of this
     * executableStep.
     *
     * @return a set of the required file structure configs, must not be null.
     */
    protected abstract Set<AbstractConfig<File>> getRequiredFileStructure();

    /**
     * Get the file configs whose value files are created during execution of this executableStep
     *
     * @return a set of the created file structure configs, must not be null.
     */
    public abstract Set<GeneratedFileStructure> getCreatedFileStructure();

    /**
     * Get the configs that are mandatory for execution of this executableStep.
     *
     * @return a set of the required configs, must not be null.
     */
    protected abstract Set<AbstractConfig<?>> getRequiredConfigs();

    /**
     * Get the configs that are not mandatory for execution of this executableStep but influence the outcome.
     * Generally if a Config.isSet() check takes place before config value usage, it is an optional config.
     *
     * @return a set of the optional configs, must not be null.
     */
    protected Set<AbstractConfig<?>> getOptionalConfigs()
    {
        return new HashSet<>();
    }

    /**
     * Update input directories if the input directory depends on configs.
     * This is necessary for the first steps of the pipeline, e.g. checkChromosomeAnnotation, MixOptions, ...
     */
    protected void updateInputDirectory()
    {
    }

    /**
     * Add all the created file structures to the TFPRIO.createdFileStructure set.
     * <p>
     * Logs an error if an already existing file structure is written to, since this would mean an in-place
     * manipulation of the file structure. This does not apply if developmentMode is active.
     */
    private void broadcastCreatedFileStructure()
    {
        for (AbstractConfig<File> createdStructure : getCreatedFileStructure())
        {
            if (TFPRIO.createdFileStructure.contains(createdStructure) && !TFPRIO.developmentMode)
            {
                logger.error("Writing to already existing structure: " + createdStructure.getName());
            } else
            {
                TFPRIO.createdFileStructure.add(createdStructure);
            }
        }
    }

    /**
     * Check if the file structure and config requirements of this executableStep are met.
     * <p>
     * The following requirements are checked:
     * <ul>
     *     <li>All the required file structures are contained by the TFPRIO.createdFileStructure set</li>
     *     <li>All the required configs are set</li>
     * </ul>
     * <p>
     * Logs a warning, if a single requirement is not met.
     *
     * @return true if all the requirements are met, otherwise false
     */
    private boolean checkRequirements()
    {
        boolean allGood = true;
        Set<AbstractConfig<?>> requiredConfigs = getRequiredConfigs();
        Set<AbstractConfig<File>> requiredFileStructure = getRequiredFileStructure();

        for (AbstractConfig<?> config : requiredConfigs)
        {
            if (!config.isSet())
            {
                allGood = false;
                logger.warn("A required config is not set: " + config.getName());
            }
        }

        for (AbstractConfig<File> fileConfig : requiredFileStructure)
        {
            if (!TFPRIO.createdFileStructure.contains(fileConfig))
            {
                allGood = false;
                if (fileConfig.get() == null)
                {
                    logger.warn("Required FileStructure is null: " + fileConfig.getName());
                } else
                {
                    String reason = ((GeneratedFileStructure) fileConfig).getNoGenerationReason();
                    logger.warn("Required FileStructure has not been created: " + fileConfig.getName() + " because: " +
                            reason);
                }
            }
        }
        return allGood;
    }

    /**
     * Delete all the file structures created by this executableStep.
     */
    private void deleteAllOutputs()
    {
        Set<GeneratedFileStructure> outputFiles = getCreatedFileStructure();
        for (GeneratedFileStructure outputFile : outputFiles)
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

    /**
     * Hash the input directories of this executableStep.
     *
     * @return hash string
     * @throws IOException if the input directories cannot be read.
     */
    private String hashInputs() throws IOException
    {
        Set<AbstractConfig<File>> requiredFileStructure = getRequiredFileStructure();
        ArrayList<File> requiredFiles = new ArrayList<>();

        for (AbstractConfig<File> fileConfig : requiredFileStructure)
        {
            requiredFiles.add(fileConfig.get());
        }

        return hashFiles(requiredFiles);
    }

    /**
     * Hash the required configs of this executableStep.
     * Required configs have to be set in order to ensure the functionality of this executableStep
     *
     * @return hash string
     */
    private String hashRequiredConfigs()
    {
        Set<AbstractConfig<?>> requiredConfigs = getRequiredConfigs();
        return util.Hashing.hashConfigs(requiredConfigs);
    }

    /**
     * Hash the optional configs of this executableStep.
     * Optional configs are configs which are checked if they are set, before their value is being used.
     *
     * @return hash string
     */
    private String hashOptionalConfigs()
    {
        Set<AbstractConfig<?>> optionalConfigs = getOptionalConfigs();
        return util.Hashing.hashConfigs(optionalConfigs);
    }

    /**
     * Hash the created file structures.
     *
     * @return hash string
     * @throws IOException if the created files cannot be read
     */
    private String hashOutputs() throws IOException
    {
        Set<GeneratedFileStructure> createdFileStructure = getCreatedFileStructure();
        ArrayList<File> createdFiles = new ArrayList<>();

        for (AbstractConfig<File> config : createdFileStructure)
        {
            createdFiles.add(config.get());
        }

        return hashFiles(createdFiles);
    }

    /**
     * Calculate the hashes of the found file structures and configs and compare them to the stored values.
     *
     * @return true if the hashes match, false if mismatch
     */
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

    /**
     * Calculate the hashes for all used configs and store them in the hash file.
     * <p>
     * Hashed config types:
     * <ul>
     *     <li>Input file structures</li>
     *     <li>Output file structures</li>
     *     <li>Required configs</li>
     *     <li>Optional configs</li>
     * </ul>
     */
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
            logger.error(e.getMessage());
        }
    }

    /**
     * Get the hash file for this executableStep.
     *
     * @return the hash file
     */
    private File getHashFile()
    {
        return extend(TFPRIO.configs.general.d_workflowHashes.get(), this.getClass().getName() + ".md5");
    }

    /**
     * Wait for all queued threads in the executorService to terminate.
     * Includes a progress and remaining time estimation.
     * TODO: Implement timout monitoring
     */
    protected void shutdown()
    {
        long total = executorService.getTaskCount();
        if (total > 0)
        {
            ExecutionTimeMeasurement timer = new ExecutionTimeMeasurement();

            executorService.shutdown();
            int timeOutMinutes = getShutDownTimeOutMinutes();

            long lastFinished = 0;
            Long latestTotalTimeGuessMillis = null;

            while (!executorService.isTerminated())
            {
                long finished = executorService.getCompletedTaskCount();
                long passedTime = timer.getDeltaMillis();

                if (finished > lastFinished)
                {
                    long millisPerStep = passedTime / finished;
                    latestTotalTimeGuessMillis = millisPerStep * total;

                    lastFinished = finished;
                }

                double percentage = 100 * (double) finished / total;
                String message = String.format("Progress: %.2f%%", percentage);

                if (lastFinished > 0)
                {
                    long remainingMillis = Math.max(latestTotalTimeGuessMillis - passedTime, 0);
                    message += ", ETA: " + ExecutionTimeMeasurement.formatMillis(remainingMillis);
                }

                logger.progress(message);
                try
                {
                    Thread.sleep(10);
                } catch (InterruptedException e)
                {
                    logger.error(e.getMessage());
                }
            }
            timer.stop();
        }
    }

    /**
     * Wait for all queued threads to terminate and reinitialize the executorService.
     */
    protected void finishAllQueuedThreads()
    {
        shutdown();
        executorService = (ThreadPoolExecutor) Executors.newFixedThreadPool(TFPRIO.configs.general.threadLimit.get());
    }

    /**
     * Iterates all class data elements in an executable step and checks if they are assigned to one of the following
     * config types:
     * <ul>
     *     <li>Created file structure</li>
     *     <li>Required file structure</li>
     *     <li>Required config</li>
     *     <li>Optional config</li>
     * </ul>
     * <p>
     * Logs a warning if a config is not assigned.
     */
    private void verifyConfigUsage()
    {
        Set<AbstractConfig<?>> registeredConfigs = new HashSet<>();
        registeredConfigs.addAll(getCreatedFileStructure());
        registeredConfigs.addAll(getRequiredConfigs());
        registeredConfigs.addAll(getOptionalConfigs());
        registeredConfigs.addAll(getRequiredFileStructure());

        for (Field field : this.getClass().getDeclaredFields())
        {
            if (AbstractConfig.class.isAssignableFrom(field.getType()))
            {
                try
                {
                    field.setAccessible(true);
                    AbstractConfig<?> current = (AbstractConfig<?>) field.get(this);
                    if (!registeredConfigs.contains(current))
                    {
                        logger.warn("Config not assigned to a category: " + field.getName() + ". May be due to a " +
                                "certain config constellation.");
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

    /**
     * Get the minutes that the executorService may take for completing all tasks.
     *
     * @return 5 if not overridden
     */
    protected int getShutDownTimeOutMinutes()
    {
        return 5;
    }

    /**
     * Get the number of threads for this ExecutableStep instance.
     *
     * @return configs->general->threadLimit if not overridden
     */
    protected int getThreadNumber()
    {
        return TFPRIO.configs.general.threadLimit.get();
    }

    /**
     * The job performed by this executableStep.
     * <p>
     * If the main job consists of multiple sub jobs that require the previous sub job to be finished, splitting the
     * process into multiple executableSteps should be considered. If this is not an option, the
     * finishAllQueuedThreads() method should be used in order to make sure that the previous sub job is finished.
     */
    protected abstract void execute();

    public void setNoExecutionReason(String noExecutionReason)
    {
        this.noExecutionReason = noExecutionReason;
    }

    public String getNoExecutionReason()
    {
        return noExecutionReason;
    }
}
