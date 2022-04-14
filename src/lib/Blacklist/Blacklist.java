package lib.Blacklist;

import tfprio.TFPRIO;
import lib.BinaryTree.ChromosomeBlacklistTrees;
import lib.ExecutableStep;
import lib.Region;
import tfprio.Workflow;
import util.Configs.ConfigTypes.AbstractConfig;
import util.FileFilters.Filters;

import java.io.*;
import java.util.*;

import static util.FileManagement.*;

public class Blacklist extends ExecutableStep
{
    private AbstractConfig<File> d_input;
    private final AbstractConfig<File> d_output = TFPRIO.configs.blacklist.fileStructure.d_newInput;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(List.of(d_input));
    }

    @Override public Set<AbstractConfig<File>> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(d_output));
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>();
    }

    @Override protected void updateInputDirectory()
    {
        d_input = Workflow.getLatestInputDirectory();
    }

    public void execute()
    {
        logger.info("Used input: " + d_input);
        logger.info("Create chromosome binary trees.");
        //CREATE BINARY TREES
        ChromosomeBlacklistTrees chromosomeTrees =
                new ChromosomeBlacklistTrees(TFPRIO.configs.blacklist.bedFilePath.get());

        logger.info("Filter input files for blacklisted regions.");
        //now filter all files

        for (File d_group : Objects.requireNonNull(d_input.get().listFiles(Filters.directoryFilter)))
        {
            File d_output_group = extend(d_output.get(), d_group.getName());

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
    }
}
