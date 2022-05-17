package util;

import tfprio.tfprio.TFPRIO;

import java.io.*;
import java.util.List;

import static util.FileManagement.extend;
import static util.FileManagement.writeFile;
import static util.ScriptExecution.executeAndWait;

/**
 * Allows generation and headless execution of igv batch files.
 */
public class IGV_Headless
{
    private final StringBuilder commandBuilder = new StringBuilder();
    private final String name;
    private final Logger logger;

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
        addCommand("exit");

        File batchFile = extend(workingDirectory, name + ".bat");

        save(batchFile);

        String command =
                "xvfb-run -s \"-terminate\" " + TFPRIO.configs.igv.pathToIGV.get().getAbsolutePath() +
                        "/igv.sh" + " -b " + batchFile.getAbsolutePath();

        executeAndWait(command, logger);
    }

    /**
     * Create a new session xml file
     *
     * @param loadFiles   a list of absolute paths to load
     * @param tdfFiles    the tdf files
     */
    public void createSession(List<String> loadFiles, List<File> tdfFiles, File f_session)
    {
        addCommand("genome " + TFPRIO.configs.igv.speciesReferenceGenome.get());

        addCommand(String.join("\nload ", loadFiles));

        //remodel tdf files if available

        for (File f_tdf : tdfFiles)
        {
            addCommand("setLogScale " + f_tdf.getName());
            addCommand("setDataRange auto " + f_tdf.getName());
            addCommand("setTrackHeight " + f_tdf.getName() + " 60");
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
        addCommand("saveSession " + f_session.getAbsolutePath());
    }
}
