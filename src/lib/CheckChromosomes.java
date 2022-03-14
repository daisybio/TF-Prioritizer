package lib;

import com2pose.COM2POSE;
import util.FileFilters.DirectoryFilter;
import util.FileFilters.FileFilter;
import util.Logger;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import static util.FileManagement.*;

public class CheckChromosomes
{
    /**
     * Make sure that chromosomes are annotated without a "chr" prefix
     */
    public void execute() throws IOException
    {
        ExecutorService executorService = Executors.newFixedThreadPool(COM2POSE.configs.general.threadLimit.get());
        Logger logger = new Logger("CHECK CHROMOSOMES");
        logger.info("Check chromosome for naming convention and alter if necessary!");

        File file_root_input = COM2POSE.configs.tepic.inputDirectory.get();
        File chr_annotation_output = COM2POSE.configs.mixOptions.fileStructure.d_preprocessingCheckChr.get();

        // I do not really know why this assignments have to take place
        //options_intern.tepic_input_prev = options_intern.tepic_input_directory;
        //options_intern.tepic_input_directory = chr_annotation_output.getAbsolutePath();

        for (File fileDir : Objects.requireNonNull(file_root_input.listFiles(new DirectoryFilter())))
        {
            String group = fileDir.getName();
            COM2POSE.groupsToHms.put(group, new HashSet<>());

            File d_outputGroup = extend(chr_annotation_output, group);

            for (File d_hm : Objects.requireNonNull(fileDir.listFiles(new DirectoryFilter())))
            {
                COM2POSE.groupsToHms.get(group).add(d_hm.getName());
                File d_outputGroupHm = extend(d_outputGroup, d_hm.getName());

                for (File f_sample : Objects.requireNonNull(d_hm.listFiles(new FileFilter())))
                {
                    executorService.execute(() ->
                    {
                        List<String> inputLines = null;
                        try
                        {
                            inputLines = readLines(f_sample);
                        } catch (IOException e)
                        {
                            logger.error("Could not read content from " + f_sample.getAbsolutePath());
                            System.exit(1);
                        }
                        StringBuilder sb_output = new StringBuilder();

                        for (String inputLine : inputLines)
                        {
                            String processedLine;
                            if (inputLine.startsWith("chr"))
                            {
                                processedLine = inputLine.substring("chr".length());
                            } else
                            {
                                processedLine = inputLine;
                            }
                            sb_output.append(processedLine);
                            sb_output.append("\n");
                        }
                        File f_outputGroupHmSample = extend(d_outputGroupHm, f_sample.getName());

                        try
                        {
                            writeFile(f_outputGroupHmSample, sb_output.toString());
                        } catch (IOException e)
                        {
                            logger.error("Could not write content to " + f_outputGroupHmSample);
                            System.exit(1);
                        }
                    });

                }
            }
        }

        executorService.shutdown();
        try
        {
            executorService.awaitTermination(5, TimeUnit.MINUTES);
        } catch (InterruptedException e)
        {
            logger.error(
                    "Chromosome annotation checking did not finish within the time limit. Message: " + e.getMessage());
            System.exit(1);
        }

        try
        {
            COM2POSE.configs.general.latestInputDirectory.setValue(chr_annotation_output);
        } catch (IllegalAccessException e)
        {
            logger.error(e.getMessage());
        }

        logger.logLine("Chromosome checking done.");
    }
}
