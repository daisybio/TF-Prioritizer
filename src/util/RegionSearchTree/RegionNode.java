package util.RegionSearchTree;

import lib.Region;

import java.util.*;

public class RegionNode
{
    private RegionNode lower = null, higher = null;
    public Region value;

    public RegionNode(Region value)
    {
        this.value = value;
    }

    public RegionNode(Iterable<Region> values)
    {
        addAllOptimized(values);
    }

    public boolean contains(Region value)
    {
        return getMatchingChild(value) != null;
    }

    public boolean matches(Region term)
    {
        return !(Math.max(term.getStart(), term.getEnd()) < Math.min(value.getStart(), value.getEnd()) ||
                Math.min(term.getStart(), term.getEnd()) > Math.max(value.getStart(), value.getEnd()));
    }

    public Region getMatchingChild(Region searchValue)
    {
        if (matches(searchValue))
        {
            return value;
        } else
        {
            int compareValue = value.compareTo(searchValue);
            if (compareValue < 0)
            {
                if (higher != null)
                {
                    return higher.getMatchingChild(searchValue);
                } else
                {
                    return null;
                }
            } else if (compareValue > 0)
            {
                if (lower != null)
                {
                    return lower.getMatchingChild(searchValue);
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

    public void add(RegionNode node)
    {
        add(node, false);
    }

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

            while (this.higher.value != null && (this.higher.value.getStart() <= this.value.getEnd()))
            {
                this.value.merge(this.higher.value);
                this.higher = this.higher.higher;
            }

            while (this.lower.value != null && (this.lower.value.getEnd() >= this.value.getStart()))
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

    public void addAllOptimized(Iterable<Region> values)
    {
        addAllOptimized(values, false);
    }

    public void addAllOptimized(Iterable<Region> values, boolean merge)
    {
        HashSet<Region> valuesSet = new HashSet<>();
        values.forEach(valuesSet::add);
        valuesSet.addAll(getAllValues());

        ArrayList<Region> valuesList = new ArrayList<>(valuesSet);

        Collections.sort(valuesList);

        addOptimizedRecursive(valuesList, merge);
    }

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

    public List<Region> getAllValuesSorted()
    {
        List<Region> list = new ArrayList<>(getAllValues());
        Collections.sort(list);
        return list;
    }

    public Set<Region> getAllValues()
    {
        Set<Region> values = new HashSet<>();
        getAllValuesRecursive(values, this);
        return values;
    }

    private void getAllValuesRecursive(Set<Region> values, RegionNode parent)
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

    public String toString()
    {
        return toStringLevel(0);
    }

    private String toStringLevel(int level)
    {
        String output = "";

        output += higher == null ? "\t".repeat(level + 1) + "NULL\n" : higher.toStringLevel(level + 1);

        output += "\t".repeat(level) + value + "\n";

        output += lower == null ? "\t".repeat(level + 1) + "NULL\n" : lower.toStringLevel(level + 1);

        return output;
    }

    public int getDepth()
    {
        int higherDepth = higher == null ? 0 : higher.getDepth();
        int lowerDepth = lower == null ? 0 : lower.getDepth();

        return Math.max(higherDepth, lowerDepth) + 1;
    }
}
