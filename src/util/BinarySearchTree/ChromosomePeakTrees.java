package util.BinarySearchTree;

import lib.Peak;
import lib.Region;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class ChromosomePeakTrees extends ChromosomeRegionTrees
{
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

    public Double getAveragePeakScore(String chromosome)
    {
        if (!hasChromosome(chromosome))
        {
            return null;
        }
        return ((PeakNode) chromosomeTrees.get(chromosome)).getAveragePeakScore();
    }

    @Override protected Node<Region> getNewNode(Iterable<Region> regions)
    {
        return new PeakNode(regions);
    }
}
