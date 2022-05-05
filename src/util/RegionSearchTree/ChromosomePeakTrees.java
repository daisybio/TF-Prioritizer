package util.RegionSearchTree;

import util.Regions.PeakRegion;
import util.Regions.Region;

import java.io.*;
import java.util.*;

/**
 * Stores a binary search tree for {@link PeakRegion} objects for each chromosome.
 */
public class ChromosomePeakTrees extends ChromosomeRegionTrees
{
    /**
     * Builds a new search tree based on a peak file
     *
     * @param source the file containing the peaks
     */
    public ChromosomePeakTrees(File source) throws IOException
    {
        addFile(source);
    }

    /**
     * Add a peak file to this object
     *
     * @param source the file to add to read the input data from
     */
    public void addFile(File source) throws IOException
    {
        Set<Region> peaks = new HashSet<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(source)))
        {
            String inputLine;
            while ((inputLine = reader.readLine()) != null)
            {
                peaks.add(new PeakRegion(inputLine));
            }
        }
        addAllOptimized(peaks);
    }

    /**
     * Get the average peak score for a given chromosome.
     *
     * @param chromosome the chromosome to get the average peak score for
     * @return the average peak score for the given chromosome
     */
    public Double getAveragePeakScore(String chromosome)
    {
        if (!hasChromosome(chromosome))
        {
            return null;
        }
        return ((PeakNode) chromosomeTrees.get(chromosome)).getAveragePeakScore();
    }

    /**
     * Create a new RegionNode based on multiple regions
     *
     * @param regions the regions to consider
     * @return the new RegionNode
     */
    @Override protected RegionNode getNewNode(Iterable<Region> regions)
    {
        return new PeakNode(regions);
    }
}
