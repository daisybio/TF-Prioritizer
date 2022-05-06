package util.RegionSearchTree;

import util.Regions.Region;

import java.util.*;

/**
 * A binary tree node that represents a certain {@link Region}.
 */
public class RegionNode
{
    /**
     * The lower child node or this node.
     */
    private RegionNode lower = null;
    /**
     * The higher child node or this node.
     */
    private RegionNode higher = null;
    /**
     * The region object represented by this node.
     */
    public Region value;

    /**
     * Creates a new instance from a single {@link Region}
     *
     * @param value the {@link Region} that is going to be represented by this node.
     */
    public RegionNode(Region value)
    {
        this.value = value;
    }

    /**
     * Creates a new instance from a list of {@link Region}
     * <p>
     * The structure of the tree is optimized in order to keep the depth as low as possible.
     *
     * @param values the regions that should be considered.
     */
    public RegionNode(Iterable<Region> values)
    {
        this(values, false);
    }

    /**
     * Creates a new instance from a list of {@link Region} with an option of region merging.
     *
     * @param values the regions to add
     * @param merge  indicates if merging should be performed
     */
    public RegionNode(Iterable<Region> values, boolean merge)
    {
        addAllOptimized(values, merge);
    }

    /**
     * Check if the node contains an overlap to the given region.
     *
     * @param value the region to search for
     * @return true if the node contains an overlap, otherwise false
     */
    public boolean hasOverlapping(Region value)
    {
        return getOverlappingChild(value) != null;
    }

    /**
     * Check if the {@link Region} represented by this node overlaps a given region.
     *
     * @param term the region to check
     * @return true if the node overlaps the given region, otherwise false
     */
    public boolean overlaps(Region term)
    {
        return value.overlaps(term);
    }

    /**
     * Get the node which overlaps a given region.
     *
     * @param searchValue the region to search for an overlap
     * @return null if no overlap is found, otherwise the found overlap
     */
    public Region getOverlappingChild(Region searchValue)
    {
        if (overlaps(searchValue))
        {
            return value;
        } else
        {
            int compareValue = value.compareTo(searchValue);
            if (compareValue < 0)
            {
                if (higher != null)
                {
                    return higher.getOverlappingChild(searchValue);
                } else
                {
                    return null;
                }
            } else if (compareValue > 0)
            {
                if (lower != null)
                {
                    return lower.getOverlappingChild(searchValue);
                } else
                {
                    return null;
                }
            } else
            {
                throw new IllegalArgumentException("Something awkward happened");
            }
        }
    }

    /**
     * Add a new node to this node.
     *
     * @param node the node to add
     */
    public void add(RegionNode node)
    {
        add(node, false);
    }

    /**
     * Add a new node to this node.
     *
     * @param node  the node to add
     * @param merge if this option is enabled and the node to add overlaps an existing one, these two nodes are merged
     */
    public void add(RegionNode node, boolean merge)
    {
        if (this.value == null)
        {
            this.value = node.value;
            return;
        }

        if (merge && this.value.overlaps(node.value))
        {
            this.value.merge(node.value);

            // Merge all higher child nodes which overlap the newly merged region
            while (this.higher != null && this.value.overlaps(this.higher.value))
            {
                this.value.merge(this.higher.value);
                this.higher = this.higher.higher;
            }

            // Merge all lower child nodes which overlap the newly merged region
            while (this.lower != null && (this.lower.value.getEnd() >= this.value.getStart()))
            {
                this.value.merge(this.lower.value);
                this.lower = this.lower.lower;
            }
            return;
        }

        if (this.value.compareTo(node.value) < 0)
        {
            if (higher != null)
            {
                higher.add(node);
            } else
            {
                higher = node;
            }
        } else
        {
            if (lower != null)
            {
                lower.add(node);
            } else
            {
                lower = node;
            }
        }
    }

    /**
     * Add multiple new {@link Region} objects to this node.
     * The structure will be optimized in order to keep the depth as low as possible.
     *
     * @param values the new regions to merge
     * @param merge  if this option is enabled, overlapping regions will be merged
     */
    public void addAllOptimized(Iterable<Region> values, boolean merge)
    {
        // Create a set of regions to add
        HashSet<Region> valuesSet = new HashSet<>();
        values.forEach(valuesSet::add);

        // Add the existing regions to the set
        valuesSet.addAll(getAllValues());

        // Create a sorted list of all regions to add
        ArrayList<Region> valuesList = new ArrayList<>(valuesSet);
        Collections.sort(valuesList);

        // Perform optimized tree creation
        addOptimizedRecursive(valuesList, merge);
    }

    /**
     * Recursively adds regions to this node by splitting the regions to add in a way that keeps the depth low.
     *
     * @param elements the regions to add
     * @param merge    indicates if overlapping regions should be merged
     */
    private void addOptimizedRecursive(List<Region> elements, boolean merge)
    {
        if (elements.isEmpty())
        {
            return;
        }

        int center = elements.size() / 2;

        List<Region> lower = elements.subList(0, center);
        List<Region> higher = elements.subList(center + 1, elements.size());

        add(new RegionNode(elements.get(center)), merge);
        addOptimizedRecursive(lower, merge);
        addOptimizedRecursive(higher, merge);
    }

    /**
     * Get a sorted list of all regions represented by this node and its sub nodes.
     *
     * @return the sorted list
     */
    public List<Region> getAllValuesSorted()
    {
        List<Region> list = new ArrayList<>(getAllValues());
        Collections.sort(list);
        return list;
    }

    /**
     * Get a set of all regions represented by this node and its sub nodes.
     *
     * @return the set of regions
     */
    public Set<Region> getAllValues()
    {
        Set<Region> values = new HashSet<>();
        getAllValuesRecursive(values, this);
        return values;
    }

    /**
     * Recursively adds all represented regions of the given parent node and its sub nodes to a given set.
     *
     * @param values the set storing the represented regions
     * @param parent the node whose represented regions should be added
     */
    private static void getAllValuesRecursive(Set<Region> values, RegionNode parent)
    {
        if (parent.value != null)
        {
            values.add(parent.value);
        }
        if (parent.lower != null)
        {
            getAllValuesRecursive(values, parent.lower);
        }
        if (parent.higher != null)
        {
            getAllValuesRecursive(values, parent.higher);
        }
    }

    @Override public String toString()
    {
        return toStringLevel(0);
    }

    /**
     * Get a formatted representation of the structure based on the given level.
     *
     * @param level the indentation level
     * @return the formatted string representing the structure
     */
    private String toStringLevel(int level)
    {
        String output = "";

        output += higher == null ? "\t".repeat(level + 1) + "NULL\n" : higher.toStringLevel(level + 1);

        output += "\t".repeat(level) + value + "\n";

        output += lower == null ? "\t".repeat(level + 1) + "NULL\n" : lower.toStringLevel(level + 1);

        return output;
    }

    /**
     * Get the depth of the structure
     *
     * @return the depth
     */
    public int getDepth()
    {
        int higherDepth = higher == null ? 0 : higher.getDepth();
        int lowerDepth = lower == null ? 0 : lower.getDepth();

        return Math.max(higherDepth, lowerDepth) + 1;
    }
}
