package org.exbio.tfprio.lib;

import org.apache.logging.log4j.Logger;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

import static org.exbio.pipejar.util.FileManagement.*;
import static org.exbio.pipejar.util.ScriptExecution.execute;
import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;


/**
 * Allows generation and headless execution of igv batch files.
 */
public class IGV_Headless {
    private static Integer xServerNum = null;
    private static Process xServer = null;
    private final StringBuilder commandBuilder = new StringBuilder();
    private final Logger logger;
    private final File igvExecutable;
    private final File igvCacheDirectory;
    private final String genome;
    private final File sessionFile;
    private final File batchFile;

    /**
     * Create a new class instance.
     *
     * @param logger the {@link org.exbio.pipejar.pipeline.ExecutableStep} logger
     */
    public IGV_Headless(Logger logger, String genome, File igvExecutable, File igvCacheDirectory,
                        File workingDirectory) {
        this.logger = logger;
        this.genome = genome;
        this.igvExecutable = igvExecutable;
        this.sessionFile = new File(workingDirectory, "session.xml");
        this.batchFile = new File(workingDirectory, "run.bat");
        this.igvCacheDirectory = igvCacheDirectory;
        addCommand("new");
        addCommand("genome " + genome);
        addCommand("load " + sessionFile.getAbsolutePath());
    }

    private static synchronized void startXServer(Logger logger) throws IOException {
        if (xServer == null) {
            xServerNum = 20;
            boolean successful = false;

            while (!successful) {
                xServerNum++;
                xServer = execute("Xvfb :" + xServerNum + " -screen 0 1131x500x16", true);
                successful = true;
                logger.info("Started XServer with ID: " + xServerNum);
            }
        }
    }

    public static synchronized void stopXServer() {
        if (xServer != null) {
            xServer.destroy();
        }
    }

    /**
     * Add a single command to an igv batch file.
     *
     * @param command the command to append, without line break.
     */
    public void addCommand(String command) {
        commandBuilder.append(command).append("\n");
    }

    private void save() {
        try {
            writeFile(batchFile, commandBuilder.toString());
        } catch (IOException e) {
            logger.error(e.getMessage());
        }
    }

    /**
     * Save and execute the created command sequence.
     */
    public void run() throws IOException {
        startXServer(logger);

        addCommand("exit");

        save();

        String command = igvExecutable.getAbsolutePath() + " -b " + batchFile.getAbsolutePath() + " --igvDirectory " +
                igvCacheDirectory.getAbsolutePath();
        
        executeAndWait(command, new HashMap<>() {{
            put("DISPLAY", ":" + xServerNum);
        }}, false);
    }

    /**
     * Create a new session xml file
     */
    public void createSession(List<File> firstPanel, List<File> secondPanel, Map<File, String> descriptions)
            throws IOException {
        if (sessionFile.exists()) {
            logger.warn("Session file already exists. Overwriting.");
            deleteFileStructure(sessionFile);
        }

        int peakHeight = 40;
        int signalHeight = 60;

        makeSureFileExists(sessionFile);

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(sessionFile))) {
            writer.write("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>");
            writer.newLine();
            writer.write("<Session genome=\"" + genome +
                    "\" hasGeneTrack=\"true\" hasSequenceTrack=\"true\" locus=\"All\" version=\"8\">");
            writer.newLine();
            writer.write("\t<Resources>");
            writer.newLine();
            Stream.of(firstPanel, secondPanel).flatMap(Collection::stream).forEach(f -> {
                try {
                    writer.write("\t\t<Resource path=\"" + f.getAbsolutePath() + "\"/>");
                    writer.newLine();
                } catch (IOException e) {
                    logger.error(e.getMessage());
                }
            });
            writer.write("\t</Resources>");
            writer.newLine();
            writer.write("\t<Panel height=\"" + ((firstPanel.size() + 1) * signalHeight) +
                    "\" name=\"DataPanel\" width=\"1131\">");
            writer.newLine();

            writer.write("\t\t<Track attributeKey=\"Reference sequence\" clazz=\"org.broad.igv.track.SequenceTrack\" " +
                    "fontSize=\"10\" id=\"Reference sequence\" name=\"Reference sequence\" sequenceTranslationStrandValue=\"POSITIVE\" " +
                    "shouldShowTranslation=\"false\" visible=\"true\"/>");
            writer.newLine();

            writer.write(
                    "\t\t<Track attributeKey=\"Refseq Genes\" clazz=\"org.broad.igv.track.FeatureTrack\" fontSize=\"10\" " +
                            "height=\"" + signalHeight +
                            "\" groupByStrand=\"false\" id=\"https://s3.amazonaws.com/igv.org.genomes/" + genome +
                            "/ncbiRefSeq.sorted.txt.gz\" " + "name=\"Refseq Genes\" visible=\"true\"/>\n");

            writer.newLine();

            firstPanel.forEach(f -> {
                try {
                    writer.write("\t\t<Track attributeKey=\"" + f.getName() +
                            "\" autoScale=\"true\" clazz=\"org.broad.igv.track.DataSourceTrack\" fontSize=\"10\" height=\"" +
                            signalHeight + "\" id=\"" + f.getAbsolutePath() + "\" name=\"" +
                            descriptions.getOrDefault(f, f.getName()) +
                            "\" renderer=\"BAR_CHART\" visible=\"true\" windowFunction=\"mean\"/>");
                    writer.newLine();
                } catch (IOException e) {
                    logger.error(e.getMessage());
                }
            });
            writer.write("\t</Panel>");
            writer.newLine();
            writer.write("\t<Panel height=\"" + (secondPanel.size() * 40) + "\" name=\"FeaturePanel\" width=\"1131\">");
            writer.newLine();
            secondPanel.forEach(f -> {
                try {
                    writer.write("\t\t<Track attributeKey=\"" + f.getName() +
                            "\" autoScale=\"true\" clazz=\"org.broad.igv.track.FeatureTrack\" fontSize=\"10\" height=\"" +
                            peakHeight + "\" id=\"" + f.getAbsolutePath() + "\" name=\"" +
                            descriptions.getOrDefault(f, f.getName()) + "\" visible=\"true\" gffTags=\"on\"/>");
                    writer.newLine();
                } catch (IOException e) {
                    logger.error(e.getMessage());
                }
            });
            writer.write("\t</Panel>");
            writer.newLine();
            writer.write("</Session>");
        }
    }
}
