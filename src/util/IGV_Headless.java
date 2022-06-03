package util;

import tfprio.tfprio.TFPRIO;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import static util.FileManagement.extend;
import static util.FileManagement.writeFile;
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
        addCommand("genome " + TFPRIO.configs.igv.speciesReferenceGenome.get());
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
     * @param loadFiles a list of files to load
     */
    public void createSession(Iterable<File> loadFiles, File f_session)
    {
        Set<String> acceptedExtensions = new HashSet<>()
        {{
            add("tdf");
            add("bed");
            add("broadPeak");
        }};
        Set<File> tdfFiles = new HashSet<>();
        Set<File> otherFiles = new HashSet<>();

        StringBuilder sb_session = new StringBuilder("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n");
        sb_session.append("<Session genome=\"").append(TFPRIO.configs.igv.speciesReferenceGenome.get())
                .append("\" hasGeneTrack=\"true\" hasSequenceTrack=\"true\" locus=\"All\" version=\"8\">\n");
        sb_session.append("\t<Resources>\n");

        for (File loadFile : loadFiles)
        {
            String extension = loadFile.getName().substring(loadFile.getName().lastIndexOf(".") + 1);

            if (acceptedExtensions.contains(extension))
            {
                sb_session.append("\t\t<Resource path=\"")
                        .append(f_session.getParentFile().toPath().relativize(loadFile.toPath())).append("\" ")
                        .append("type=\"").append(extension.toLowerCase()).append("\"/>\n");

                if ("tdf".equals(extension))
                {
                    tdfFiles.add(loadFile);
                } else
                {
                    otherFiles.add(loadFile);
                }
            } else
            {
                logger.warn("Unknown file type to load for igv: " + loadFile.getAbsolutePath());
            }
        }
        sb_session.append("\t</Resources>\n");

        int tdfTrackHeight = 60;
        sb_session.append("\t<Panel height=\"").append((tdfFiles.size() + 1) * tdfTrackHeight)
                .append("\" name=\"DataPanel\" width=\"1131\">\n");
        sb_session.append("\t\t<Track attributeKey=\"Reference sequence\" clazz=\"org.broad.igv.track.SequenceTrack\"" +
                " fontSize=\"10\" id=\"Reference sequence\" name=\"Reference sequence\" sequenceTranslationStrandValue=\"POSITIVE\" shouldShowTranslation=\"false\" visible=\"true\"/>\n");
        sb_session.append("\t\t<Track attributeKey=\"Refseq Genes\" clazz=\"org.broad.igv.track.FeatureTrack\" " +
                        "fontSize=\"10\" height=\"" + tdfTrackHeight + "\" groupByStrand=\"false\" id=\"https://s3" +
                        ".amazonaws" + ".com/igv" + ".org.genomes/").append(TFPRIO.configs.igv.speciesReferenceGenome.get())
                .append("/ncbiRefSeq.sorted.txt.gz\" name=\"Refseq Genes\" visible=\"true\"/>\n");

        for (File file : tdfFiles)
        {
            sb_session.append("\t\t<Track attributeKey=\"").append(file.getName())
                    .append("\" autoScale=\"true\" clazz=\"org")
                    .append(".broad.igv.track.DataSourceTrack\" fontSize=\"10\" height=\"" + tdfTrackHeight + "\" " +
                            "id=\"").append(file.getAbsolutePath()).append("\" name=\"").append(file.getName())
                    .append("\" renderer=\"BAR_CHART\" ").append("visible=\"true\" windowFunction=\"mean\"/>\n");
        }

        sb_session.append("\t</Panel>\n");
        sb_session.append("\t<Panel height=\"").append(otherFiles.size() * 40)
                .append("\" name=\"FeaturePanel\" width=\"1131\">\n");

        Map<String, Set<File>> symbol_files = new HashMap<>();
        Set<String> symbolsWithPrediction = new HashSet<>();

        for (File file : otherFiles)
        {
            String fileName = file.getName().toUpperCase().substring(0, file.getName().lastIndexOf("."));
            ArrayList<String> nameParts = new ArrayList<>(List.of(fileName.split("_")));
            nameParts.removeAll(TFPRIO.existingHms.stream().map(String::toUpperCase).collect(Collectors.toList()));
            nameParts.removeAll(
                    TFPRIO.groupsToHms.keySet().stream().map(String::toUpperCase).collect(Collectors.toList()));
            nameParts.removeAll(TFPRIO.groupCombinationsToHms.keySet().stream().map(String::toUpperCase)
                    .collect(Collectors.toList()));

            if (nameParts.size() != 1)
            {
                logger.warn(String.valueOf(nameParts));
            }
            String name = nameParts.get(0);
            if (!symbol_files.containsKey(name))
            {
                symbol_files.put(name, new HashSet<>());
            }
            if (file.getAbsolutePath().startsWith(
                    TFPRIO.configs.tepic.fileStructure.d_postprocessing_trapPredictedBeds.get().getAbsolutePath()))
            {
                symbolsWithPrediction.add(name);
            }

            symbol_files.get(name).add(file);
        }

        Consumer<Set<File>> addFiles = (set) ->
        {
            Set<String> addedTrackNames = new HashSet<>();

            for (File file : set)
            {
                String trackName = file.getName();

                try (BufferedReader reader = new BufferedReader(new FileReader(file)))
                {
                    String firstLine = reader.readLine();
                    if (firstLine.startsWith("track name="))
                    {
                        trackName = firstLine.substring(firstLine.indexOf("\"") + 1, firstLine.lastIndexOf("\""));
                    }
                } catch (IOException ignore)
                {
                }

                if (addedTrackNames.contains(trackName))
                {
                    continue;
                }

                sb_session.append("\t\t<Track attributeKey=\"").append(file.getName())
                        .append("\" autoScale=\"true\" clazz=\"org")
                        .append(".broad.igv.track.FeatureTrack\" fontSize=\"10\" height=\"40\" id=\"")
                        .append(file.getAbsolutePath()).append("\" name=\"").append(trackName)
                        .append("\" visible=\"true\"/>\n");

                addedTrackNames.add(trackName);
            }
        };

        for (String symbolWithPrediction : symbolsWithPrediction.stream().sorted().collect(Collectors.toList()))
        {
            addFiles.accept(symbol_files.get(symbolWithPrediction));
        }

        for (String symbolWithoutPrediction : symbol_files.keySet().stream()
                .filter(entry -> !symbolsWithPrediction.contains(entry)).sorted().collect(Collectors.toList()))
        {
            addFiles.accept(symbol_files.get(symbolWithoutPrediction));
        }


        sb_session.append("\t</Panel>\n");
        sb_session.append("</Session>");

        writeFile(f_session, sb_session.toString(), logger);
        addCommand("load " + f_session.getAbsolutePath());

        for (File file : tdfFiles)
        {
            addCommand("setDataRange auto " + file.getName());
        }

        //include enhancer regions of interest if available
        if (TFPRIO.configs.igv.enhancerDatabases.isSet())
        {
            File f_database =
                    extend(TFPRIO.configs.deSeq2.fileStructure.f_preprocessing_genePositions_mergedEnhancerDbs.get());

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

        addCommand("saveSession " + extend(f_session.getParentFile(), "processed.xml").getAbsolutePath());
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
                xServer = execute("Xvfb :" + xServerNum + " -screen 0 1131x500x16", logger);
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
