package lib.MixOptions;

import com2pose.COM2POSE;
import lib.BinaryTree.ChromosomePeakTrees;
import lib.Peak;
import lib.Region;
import util.Comparators.ChromosomeComparator;
import util.Logger;

import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import static util.FileManagement.*;

public class MixMutuallyExclusive
{
    private final Logger logger;

    public MixMutuallyExclusive()
    {
        this.logger = new Logger("MUTUALLY-EXCLUSIVE-PEAKS");
    }

    public void execute()
    {
        ExecutorService executorService = Executors.newFixedThreadPool(COM2POSE.configs.general.threadLimit.get());
        logger.info("Start mutually exclusive peaks calculation.");

        logger.info("Preprocessing mutually exclusive peaks for binary tree comparison.");

        File d_input = COM2POSE.configs.general.latestInputDirectory.get();

        logger.info("Used data: " + d_input);

        logger.info("Starting binary tree filtering.");

        for (String group1 : COM2POSE.groupsToHms.keySet())
        {
            for (String group2 : COM2POSE.groupsToHms.keySet())
            {
                if (group1.compareTo(group2) >= 0)
                {
                    continue;
                }
                String combined = group1 + "_" + group2;

                for (String hm : COM2POSE.groupsToHms.get(group1))
                {
                    if (!COM2POSE.groupsToHms.get(group2).contains(hm))
                    {
                        continue;
                    }

                    File d_combinedOutput =
                            extend(COM2POSE.configs.mixOptions.fileStructure.d_mutuallyExclusive_input.get(), combined,
                                    hm);

                    //input files
                    File d_unfilteredInput1 = extend(d_input, group1, hm);
                    File d_unfilteredInput2 = extend(d_input, group2, hm);
                    File[] inputFiles1 = Objects.requireNonNull(d_unfilteredInput1.listFiles());
                    File[] inputFiles2 = Objects.requireNonNull(d_unfilteredInput2.listFiles());
                    assert inputFiles1.length == 1;
                    assert inputFiles2.length == 1;
                    File f_unfilteredInput1 = inputFiles1[0];
                    File f_unfilteredInput2 = inputFiles2[0];

                    //filter tp2
                    {
                        File f_output2 = extend(d_combinedOutput, f_unfilteredInput2.getName());
                        executorService.execute(() -> filter(f_unfilteredInput2, f_unfilteredInput1, f_output2));
                    }

                    //filter tp1
                    {
                        File f_output1 = extend(d_combinedOutput, f_unfilteredInput1.getName());
                        executorService.execute(() -> filter(f_unfilteredInput1, f_unfilteredInput2, f_output1));
                    }
                }
            }
        }

        executorService.shutdown();
        try
        {
            executorService.awaitTermination(5, TimeUnit.MINUTES);
        } catch (InterruptedException e)
        {
            logger.error("Binary tree filtering did not finish in time.");
            System.exit(1);
        }
        logger.info("Finished binary tree filtering.");
    }

    private void filter(File f_source, File f_filter, File f_output)
    {
        try
        {
            ChromosomePeakTrees chromosomeTreesGroup1 = new ChromosomePeakTrees(f_filter);
            makeSureFileExists(f_output);
            Set<Peak> peaks = new HashSet<>();
            try (BufferedReader reader = new BufferedReader(new FileReader(f_source)))
            {
                String inputLine;


                while ((inputLine = reader.readLine()) != null)
                {
                    peaks.add(new Peak(inputLine));
                }
            }
            Map<String, Set<Peak>> filteredPeaks = filterPeaksByChromosomeTree(peaks, chromosomeTreesGroup1);

            try (BufferedWriter writer = new BufferedWriter(new FileWriter(f_output)))
            {
                List<String> chromosomes = new ArrayList<>(filteredPeaks.keySet());
                chromosomes.sort(new ChromosomeComparator());

                for (String chromosome : chromosomes)
                {
                    ArrayList<Peak> peakList = new ArrayList<>(filteredPeaks.get(chromosome));
                    Collections.sort(peakList);
                    for (Peak peak : peakList)
                    {
                        writer.write(peak.toString());
                        writer.newLine();
                    }
                }
            }
        } catch (IOException e)
        {
            e.printStackTrace();
        }
    }

    private Map<String, Set<Peak>> filterPeaksByChromosomeTree(Iterable<Peak> peaks, ChromosomePeakTrees peakTrees)
    {
        Map<String, Set<Peak>> output = new HashMap<>();

        for (String chromosome : peakTrees.getChromosomes())
        {
            output.put(chromosome, new HashSet<>());
        }

        for (Peak peak : peaks)
        {
            if (!peakTrees.hasChromosome(peak.getChromosome()))
            {
                continue;
            }

            if (COM2POSE.configs.mixOptions.mutuallyExclusiveDifferentialPeakSignals.get())
            {
                Peak match = (Peak) peakTrees.getMatchingChild(peak);

                if (match == null)
                {
                    //this peak is mutually exclusive
                    output.get(peak.getChromosome()).add(peak);
                    continue;
                }

                double peak_difference = match.getScore() - peak.getScore();

                peak_difference = Math.abs(peak_difference);

                if (peak_difference > peakTrees.getAveragePeakScore(peak.getChromosome()))
                {
                    //this peak is mutually exclusive in peak score
                    output.get(peak.getChromosome()).add(peak);
                }
            } else
            {
                Region match = peakTrees.getMatchingChild(peak);
                if (match == null)
                {
                    //this peak is mutually exclusive
                    output.get(peak.getChromosome()).add(peak);
                }
            }
        }
        return output;
    }
}
