package util.RegionSearchTree;

import util.Regions.Region;

import java.util.*;

/**
 * Stores one binary search tree for each chromosome
 */
public class ChromosomeRegionTrees
{
    /**
     * Maps the chromosome name to its corresponding binary tree
     */
    final Map<String, RegionNode> chromosomeTrees = new HashMap<>();

    /**
     * Get an overlapping region if the chromosome tree contains such a region
     *
     * @param region the region to search for overlaps with
     * @return null if no overlap is found, otherwise the overlapping region
     */
    public Region getOverlappingChild(Region region)
    {
        if (!hasChromosome(region.getChromosome()))
        {
            return null;
        }
        return chromosomeTrees.get(region.getChromosome()).getOverlappingChild(region);
    }

    /**
     * Check if the chromosomeTrees contain an overlap for a certain region
     *
     * @param region the region to search for overlaps with
     * @return true if an overlap exists, otherwise false
     */
    public boolean hasOverlap(Region region)
    {
        return getOverlappingChild(region) != null;
    }

    /**
     * Create a new Node based on multiple regions. Can be overridden in order to make sure that the derived class
     * stays consistent.
     *
     * @param regions the regions to consider
     * @return the newly created RegionNode object
     */
    protected RegionNode getNewNode(Iterable<Region> regions)
    {
        return new RegionNode(regions);
    }

    /**
     * @return all the region objects represented by this node and its sub nodes in order
     */
    public Map<String, List<Region>> getAllRegionsSorted()
    {
        Map<String, List<Region>> chromosomeLists = new HashMap<>();

        for (String chromosome : chromosomeTrees.keySet())
        {
            chromosomeLists.put(chromosome, chromosomeTrees.get(chromosome).getAllValuesSorted());
        }
        return chromosomeLists;
    }

    /**
     * Add multiple regions to this trees while keeping the tree depth as low as possible
     *
     * @param regions the regions to add
     */
    public void addAllOptimized(Iterable<Region> regions)
    {
        addAllOptimized(regions, false);
    }

    /**
     * Add multiple regions to this trees while keeping the tree depth as low as possible.
     *
     * @param regions the regions to add
     * @param merge   indicates if overlapping regions should be merged
     */
    public void addAllOptimized(Iterable<Region> regions, boolean merge)
    {
        Map<String, Set<Region>> chromosomeSets = new HashMap<>();

        for (Region peak : regions)
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
                chromosomeTrees.put(chromosome, getNewNode(chromosomeSets.get(chromosome)));
            } else
            {
                chromosomeTrees.get(chromosome).addAllOptimized(chromosomeSets.get(chromosome), merge);
            }
        }
    }

    /**
     * @return a set of the chromosomes available in this object
     */
    public Set<String> getChromosomes()
    {
        return chromosomeTrees.keySet();
    }

    /**
     * Check if this object contains data for a given chromosome
     *
     * @param chromosome the chromosome to check for existence
     * @return true if the chromosome exists, otherwise false
     */
    public boolean hasChromosome(String chromosome)
    {
        return chromosomeTrees.containsKey(chromosome);
    }
}
