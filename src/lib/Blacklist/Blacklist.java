package lib.Blacklist;

import com2pose.COM2POSE;
import lib.BinaryTree.ChromosomeBlacklistTrees;
import lib.BinaryTree.RegionNode;
import lib.ExecutableStep;
import lib.Region;
import util.FileFilters.Filters;
import util.Logger;

import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import static util.FileManagement.*;

public class Blacklist extends ExecutableStep
{
    public void execute()
    {
        File d_input = COM2POSE.configs.general.latestInputDirectory.get();

        File d_output = COM2POSE.configs.blacklist.fileStructure.d_newInput.get();

        logger.info("Used input: " + d_input);
        logger.info("Create chromosome binary trees.");
        //CREATE BINARY TREES
        ChromosomeBlacklistTrees chromosomeTrees =
                new ChromosomeBlacklistTrees(COM2POSE.configs.blacklist.bedFilePath.get());

        logger.info("Filter input files for blacklisted regions.");
        //now filter all files

        for (File d_group : Objects.requireNonNull(d_input.listFiles(Filters.directoryFilter)))
        {
            File d_output_group = extend(d_output, d_group.getName());

            for (File d_hm : Objects.requireNonNull(d_group.listFiles(Filters.directoryFilter)))
            {
                File d_output_groupHm = extend(d_output_group, d_hm.getName());

                for (File f_sample : Objects.requireNonNull(d_hm.listFiles(Filters.fileFilter)))
                {
                    executorService.execute(() ->
                    {
                        logger.info(
                                "Filter: " + d_group.getName() + ": " + d_hm.getName() + " - " + f_sample.getName());

                        File f_output = extend(d_output_groupHm, f_sample.getName());
                        try
                        {
                            makeSureFileExists(f_output);
                            try (BufferedWriter writer = new BufferedWriter(new FileWriter(f_output));
                                 BufferedReader reader = new BufferedReader(new FileReader(f_sample)))
                            {
                                String inputLine;

                                while ((inputLine = reader.readLine()) != null)
                                {
                                    String[] split = inputLine.split("\t");

                                    String chromosome = split[0];

                                    Region region = new Region(split);

                                    if (chromosomeTrees.hasChromosome(chromosome))
                                    {
                                        if (!chromosomeTrees.contains(region))
                                        {
                                            writer.write(inputLine);
                                            writer.newLine();
                                        }
                                    }
                                }
                            }
                        } catch (IOException e)
                        {
                            logger.error(e.getMessage());
                            System.exit(1);
                        }
                    });
                }
            }
        }
        try
        {
            COM2POSE.configs.general.latestInputDirectory.setValue(d_output);
        } catch (IllegalAccessException e)
        {
            logger.error(e.getMessage());
        }
    }
}
