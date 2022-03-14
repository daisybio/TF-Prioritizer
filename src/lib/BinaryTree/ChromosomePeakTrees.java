package lib.BinaryTree;

import lib.Peak;
import lib.Region;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class ChromosomePeakTrees
{
    final Map<String, PeakNode> chromosomeTrees = new HashMap<>();

    public ChromosomePeakTrees(File source)
    {
        addFile(source);
    }

    public void addFile(File source)
    {
        Set<Region> peaks = new HashSet<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(source)))
        {
            String inputLine;
            while ((inputLine = reader.readLine()) != null)
            {
                peaks.add(new Peak(inputLine));
            }
        } catch (IOException e)
        {
            e.printStackTrace();
        }
        addAllOptimized(peaks);
    }

    public void addAllOptimized(Iterable<Region> peaks)
    {
        Map<String, Set<Region>> chromosomeSets = new HashMap<>();

        for (Region peak : peaks)
        {
            if (!chromosomeSets.containsKey(peak.getChromosome()))
            {
                chromosomeSets.put(peak.getChromosome(), new HashSet<>());
            }
            chromosomeSets.get(peak.getChromosome()).add(peak);
        }

        for (String chromosome : chromosomeSets.keySet())
        {
            if (!chromosomeTrees.containsKey(chromosome))
            {
                chromosomeTrees.put(chromosome, new PeakNode(chromosomeSets.get(chromosome)));
            } else
            {
                chromosomeTrees.get(chromosome).addAllOptimized(chromosomeSets.get(chromosome));
            }
        }
    }

    public Region getMatchingChild(Peak peak)
    {
        if (!hasChromosome(peak.getChromosome()))
        {
            return null;
        }
        return chromosomeTrees.get(peak.getChromosome()).getMatchingChild(peak);
    }

    public Map<String, List<Region>> getAllPeaksSorted()
    {
        Map<String, List<Region>> chromosomeLists = new HashMap<>();

        for (String chromosome : chromosomeTrees.keySet())
        {
            chromosomeLists.put(chromosome, chromosomeTrees.get(chromosome).getAllValuesSorted());
        }
        return chromosomeLists;
    }

    public boolean hasChromosome(String chromosome)
    {
        return chromosomeTrees.containsKey(chromosome);
    }

    public Double getAveragePeakScore(String chromosome)
    {
        if (!hasChromosome(chromosome))
        {
            return null;
        }
        return chromosomeTrees.get(chromosome).getAveragePeakScore();
    }

    public Set<String> getChromosomes()
    {
        return chromosomeTrees.keySet();
    }
}
