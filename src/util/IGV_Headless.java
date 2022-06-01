package util;

import tfprio.tfprio.TFPRIO;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.TimeUnit;

import static util.FileManagement.*;
import static util.ScriptExecution.execute;

/**
 * Allows generation and headless execution of igv batch files.
 */
public class IGV_Headless
{
    private final StringBuilder commandBuilder = new StringBuilder();
    private final String name;
    private final Logger logger;

    private static Integer xServerNum = null;
    private static Process xServer = null;

    /**
     * Create a new class instance.
     *
     * @param name   the batch file name
     * @param logger the {@link lib.ExecutableStep} logger
     */
    public IGV_Headless(String name, Logger logger)
    {
        this.logger = logger;
        this.name = name;
        addCommand("new");
    }

    /**
     * Add a single command to an igv batch file.
     *
     * @param command the command to append, without line break.
     */
    public void addCommand(String command)
    {
        commandBuilder.append(command).append("\n");
    }

    private void save(File file)
    {
        try
        {
            writeFile(file, commandBuilder.toString());
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }
    }

    /**
     * Save and execute the created command sequence.
     *
     * @param workingDirectory the directory to save the batch file in. The generated screenshots will also be saved
     *                         to this directory, if not referenced by absolute paths.
     */
    public void run(File workingDirectory)
    {
        startXServer(logger);

        addCommand("exit");

        File batchFile = extend(workingDirectory, name + ".bat");

        save(batchFile);

        String command =
                TFPRIO.configs.igv.pathToIGV.get().getAbsolutePath() + "/igv.sh -b " + batchFile.getAbsolutePath();

        boolean worked = false;
        int attempt = 0;

        while (!worked)
        {
            try
            {
                Process igv = execute(command, logger, new HashMap<>()
                {{
                    put("DISPLAY", ":" + xServerNum);
                }}, true);

                boolean returnValue = igv.waitFor(5, TimeUnit.MINUTES);
                if (returnValue)
                {
                    worked = true;
                } else
                {
                    logger.warn("IGV did not work on attempt #" + (attempt + 1));
                    attempt++;
                    igv.destroy();
                    logger.debug("Destroyed igv process.");
                    if (attempt >= 10)
                    {
                        logger.error("Exiting since igv did not work.");
                    }
                }
            } catch (InterruptedException e)
            {
                logger.error(e.getMessage());
            }
        }
    }

    /**
     * Create a new session xml file
     *
     * @param loadFiles a list of absolute paths to load
     * @param tdfFiles  the tdf files
     */
    public void createSession(List<String> loadFiles, List<File> tdfFiles, File f_session)
    {
        addCommand("maxPanelHeight 3000");
        addCommand("genome " + TFPRIO.configs.igv.speciesReferenceGenome.get());
        Set<String> nonTdfFiles = new HashSet<>();

        for (String loadFile : loadFiles)
        {
            if (loadFile.endsWith(".tdf"))
            {
                addCommand("load " + loadFile);
                File f_tdf = new File(loadFile);
                addCommand("setLogScale true " + f_tdf.getName());
                addCommand("setDataRange auto " + f_tdf.getName());
                addCommand("setTrackHeight " + f_tdf.getName() + " 60");
            } else
            {
                nonTdfFiles.add(loadFile);
            }
        }

        for (String entry : nonTdfFiles)
        {
            addCommand("load " + entry);
        }

        //include enhancer regions of interest if available
        if (TFPRIO.configs.igv.enhancerDatabases.isSet())
        {
            for (String enhancerDB : TFPRIO.configs.igv.enhancerDatabases.get())
            {
                File f_database =
                        extend(TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_genePositions_enhancerDBs.get(),
                                enhancerDB + ".bed");

                try (BufferedReader br_mergedDB = new BufferedReader(new FileReader(f_database)))
                {
                    String line_mergedDB;
                    br_mergedDB.readLine();
                    while ((line_mergedDB = br_mergedDB.readLine()) != null)
                    {
                        String[] split = line_mergedDB.split("\t");

                        addCommand("region " + split[0] + " " + split[1] + " " + split[2]);
                    }
                } catch (IOException e)
                {
                    logger.error(e.getMessage());
                }
            }
        }
        makeSureFileExists(f_session, logger);
        addCommand("saveSession " + f_session.getAbsolutePath());
    }

    private static synchronized void startXServer(Logger logger)
    {
        if (xServer == null)
        {
            xServerNum = 20;
            boolean successful = false;

            while (!successful)
            {
                xServerNum++;
                xServer = execute("Xvfb :" + xServerNum + " -screen 0 4000x4000x16", logger);
                successful = true;
                logger.info("Started XServer with ID: " + xServerNum);
            }
        }
    }

    public static void stopXServer()
    {
        if (xServer != null)
        {
            xServer.destroy();
        }
    }
}
