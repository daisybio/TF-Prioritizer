package lib.Blacklist;

import com2pose.COM2POSE;
import lib.Region;
import util.FileFilters.Filters;
import util.Logger;

import java.io.*;
import java.util.*;

import static util.FileManagement.*;

public class Blacklist
{
    private final Logger logger;

    public Blacklist()
    {
        logger = new Logger("BLACKLIST");
    }

    public void execute() throws IOException
    {
        preprocess();
        filter();
    }

    private void preprocess() throws IOException
    {
        logger.info("Start preprocessing blacklist");
        String header = "CHR\tLEFT_BORDER\tRIGHT_BORDER\tSIGNAL";

        File f_blacklist_pre_chr = COM2POSE.configs.blacklist.fileStructure.d_preprocessing_perChr.get();

        {
            ArrayList<BlacklistedRegion> chromosome_ius = new ArrayList<>();

            List<String> inputLines = readLines(COM2POSE.configs.blacklist.bedFilePath.get());

            String currentChromosome =
                    inputLines.get(0).substring(0, inputLines.get(0).indexOf("\t")).replace("chr", "");

            for (int i = 0; i < inputLines.size(); i++)
            {
                String inputLine = inputLines.get(i);
                boolean lastLine = (i == inputLines.size() - 1);
                String[] split = inputLine.split("\t");

                String chromosome = split[0].replace("chr", "");

                if (lastLine)
                {
                    BlacklistedRegion iu = new BlacklistedRegion(split);

                    if (COM2POSE.configs.blacklist.signalsToIgnore.get().contains(iu.signal))
                    {
                        chromosome_ius.add(iu);
                    }
                }

                if (!chromosome.equals(currentChromosome) || lastLine)
                {
                    File targetFile = extend(f_blacklist_pre_chr, currentChromosome + ".tsv");
                    makeSureFileExists(targetFile);

                    try (BufferedWriter chromosomeWriter = new BufferedWriter(new FileWriter(targetFile)))
                    {
                        chromosomeWriter.write(header);
                        chromosomeWriter.newLine();

                        Collections.sort(chromosome_ius);


                        for (BlacklistedRegion iu : chromosome_ius)
                        {
                            chromosomeWriter.write(iu.toString());
                            chromosomeWriter.newLine();
                        }
                    }

                    File targetFileSorted =
                            extend(COM2POSE.configs.blacklist.fileStructure.d_preprocessing_sorted.get(),
                                    currentChromosome + ".tsv");
                    makeSureFileExists(targetFileSorted);

                    try (BufferedWriter writerSort = new BufferedWriter(new FileWriter(targetFileSorted)))
                    {
                        ArrayList<BlacklistedRegion> newly_ordered = new ArrayList<>();

                        recursive_split_BL(chromosome_ius, newly_ordered);

                        writerSort.write(header);
                        writerSort.newLine();

                        for (BlacklistedRegion entry : newly_ordered)
                        {
                            writerSort.write(entry.toString());
                            writerSort.newLine();
                        }
                    }

                    chromosome_ius.clear();
                    currentChromosome = split[0].replace("chr", "");
                }

                if (!lastLine)
                {
                    BlacklistedRegion iu = new BlacklistedRegion(inputLine.split("\t"));

                    if (COM2POSE.configs.blacklist.signalsToIgnore.get().contains(iu.signal))
                    {
                        chromosome_ius.add(iu);
                    }
                }
            }
        }
        logger.info("Finished preprocessing blacklist");
    }

    private void filter()
    {
        File d_input;

        if (COM2POSE.configs.tepic.tfBindingSiteSearch.get().equals("BETWEEN") ||
                COM2POSE.configs.tepic.tfBindingSiteSearch.get().equals("EXCL_BETWEEN"))
        {
            d_input = COM2POSE.configs.mixOptions.fileStructure.d_footprintsBetweenPeaks.get();
        } else if (COM2POSE.configs.mixOptions.level.get().equals("SAMPLE_LEVEL"))
        {
            d_input = COM2POSE.configs.mixOptions.fileStructure.d_sampleMix.get();
        } else
        {
            d_input = COM2POSE.configs.mixOptions.fileStructure.d_preprocessingCheckChr.get();
        }
        File d_output = COM2POSE.configs.blacklist.fileStructure.d_newInput.get();

        logger.info("Used input: " + d_input);
        logger.info("Create chromosome binary trees.");
        //CREATE BINARY TREES
        HashMap<String, BinaryTreeNode> chromosomeTree = new HashMap<>();

        File d_blacklistInput = COM2POSE.configs.blacklist.fileStructure.d_preprocessing_sorted.get();

        for (File f_chromosome : Objects.requireNonNull(d_blacklistInput.listFiles(Filters.fileFilter)))
        {
            String chromosome = f_chromosome.getName().substring(0, f_chromosome.getName().lastIndexOf("."));

            ArrayList<BlacklistedRegion> region = new ArrayList<>();

            try (BufferedReader blacklistReader = new BufferedReader(new FileReader(f_chromosome)))
            {
                String line_chr;
                blacklistReader.readLine();

                while ((line_chr = blacklistReader.readLine()) != null)
                {
                    BlacklistedRegion iu = new BlacklistedRegion(line_chr.split("\t"));
                    region.add(iu);
                }
            } catch (IOException e)
            {
                e.printStackTrace();
            }

            BinaryTreeNode tree = new BinaryTreeNode(region.get(0));

            for (int i = 1; i < region.size(); i++)
            {
                tree.add(region.get(i));
            }

            chromosomeTree.put(chromosome, tree);
        }

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
                    logger.info("Filter: " + d_group.getName() + ": " + d_hm.getName() + " - " + f_sample.getName());
                    //FILTER HERE WITH BINARY TREE

                    File f_output = extend(d_output_groupHm, f_sample.getName());
                    try
                    {
                        makeSureFileExists(f_output);
                    } catch (IOException e)
                    {
                        logger.error(e.getMessage());
                        System.exit(1);
                    }

                    try (BufferedWriter writer = new BufferedWriter(new FileWriter(f_output));
                         BufferedReader reader = new BufferedReader(new FileReader(f_sample)))
                    {
                        String inputLine;

                        while ((inputLine = reader.readLine()) != null)
                        {
                            String[] split = inputLine.split("\t");

                            String chromosome = split[0];

                            Region region =
                                    new Region(chromosome, Integer.parseInt(split[1]), Integer.parseInt(split[2]));

                            if (chromosomeTree.containsKey(chromosome))
                            {
                                BinaryTreeNode tree = chromosomeTree.get(chromosome);

                                if (!tree.contains(region))
                                {
                                    writer.write(inputLine);
                                    writer.newLine();
                                }
                            }
                        }
                    } catch (IOException e)
                    {
                        e.printStackTrace();
                    }
                }
            }
        }
        logger.info("Finished filtering for blacklisted regions.");
    }

    private void recursive_split_BL(ArrayList<BlacklistedRegion> region, ArrayList<BlacklistedRegion> newly_ordered)
    {
        if (region.size() > 0)
        {
            BlacklistedRegion median = region.get(region.size() / 2);
            newly_ordered.add(median);

            ArrayList<BlacklistedRegion> region_left = new ArrayList<>();

            for (int i = 0; i < region.size() / 2; i++)
            {
                region_left.add(region.get(i));
            }
            ArrayList<BlacklistedRegion> region_right = new ArrayList<>();
            for (int i = region.size() / 2 + 1; i < region.size(); i++)
            {
                region_right.add(region.get(i));
            }

            recursive_split_BL(region_left, newly_ordered);
            recursive_split_BL(region_right, newly_ordered);
        }
    }
}
