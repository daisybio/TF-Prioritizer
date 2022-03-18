package lib.MixOptions;

import tfprio.TFPRIO;
import lib.BinaryTree.ChromosomePeakTrees;
import lib.ExecutableStep;
import lib.Peak;
import lib.Region;
import util.Comparators.ChromosomeComparator;
import util.Configs.Config;

import java.io.*;
import java.util.*;

import static util.FileManagement.*;

public class MixMutuallyExclusive extends ExecutableStep
{
    private final Config<File> d_input = TFPRIO.latestInputDirectory;
    private final Config<File> d_output = TFPRIO.configs.mixOptions.fileStructure.d_mutuallyExclusive_input;
    private final Config<Boolean> mutuallyExclusiveDifferentialPeakSignals =
            TFPRIO.configs.mixOptions.mutuallyExclusiveDifferentialPeakSignals;

    @Override protected Set<Config<File>> getRequiredFileStructure()
    {
        return new HashSet<>(List.of(d_input));
    }

    @Override protected Set<Config<File>> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(d_output));
    }

    @Override protected Set<Config<?>> getRequiredConfigs()
    {
        return new HashSet<>(List.of(mutuallyExclusiveDifferentialPeakSignals));
    }

    public void execute()
    {
        logger.info("Start mutually exclusive peaks calculation.");

        logger.info("Preprocessing mutually exclusive peaks for binary tree comparison.");

        logger.info("Used data: " + d_input);

        logger.info("Starting binary tree filtering.");

        for (String group1 : TFPRIO.groupsToHms.keySet())
        {
            for (String group2 : TFPRIO.groupsToHms.keySet())
            {
                if (group1.compareTo(group2) >= 0)
                {
                    continue;
                }
                String combined = group1 + "_" + group2;

                for (String hm : TFPRIO.groupsToHms.get(group1))
                {
                    if (!TFPRIO.groupsToHms.get(group2).contains(hm))
                    {
                        continue;
                    }

                    File d_combinedOutput = extend(d_output.get(), combined, hm);

                    //input files
                    File d_unfilteredInput1 = extend(d_input.get(), group1, hm);
                    File d_unfilteredInput2 = extend(d_input.get(), group2, hm);
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

            if (mutuallyExclusiveDifferentialPeakSignals.get())
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
