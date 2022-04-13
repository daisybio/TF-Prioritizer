package lib.BinaryTree;

import lib.Region;
import util.BinarySearchTree.Node;

import java.util.*;

public class ChromosomeRegionTrees
{
    final Map<String, Node<Region>> chromosomeTrees = new HashMap<>();

    public Region getMatchingChild(Region region)
    {
        if (!hasChromosome(region.getChromosome()))
        {
            return null;
        }
        return chromosomeTrees.get(region.getChromosome()).getMatchingChild(region);
    }

    public boolean contains(Region region)
    {
        return getMatchingChild(region) != null;
    }

    protected Node<Region> getNewNode(Region region)
    {
        return new RegionNode(region);
    }

    protected Node<Region> getNewNode(Iterable<Region> regions)
    {
        return new RegionNode(regions);
    }

    public Map<String, List<Region>> getAllRegionsSorted()
    {
        Map<String, List<Region>> chromosomeLists = new HashMap<>();

        for (String chromosome : chromosomeTrees.keySet())
        {
            chromosomeLists.put(chromosome, chromosomeTrees.get(chromosome).getAllValuesSorted());
        }
        return chromosomeLists;
    }

    public void addAllOptimized(Iterable<Region> regions)
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
                chromosomeTrees.get(chromosome).addAllOptimized(chromosomeSets.get(chromosome));
            }
        }
    }

    public Set<String> getChromosomes()
    {
        return chromosomeTrees.keySet();
    }

    public boolean hasChromosome(String chromosome)
    {
        return chromosomeTrees.containsKey(chromosome);
    }
}
