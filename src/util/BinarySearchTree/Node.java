package util.BinarySearchTree;

import lib.Region;

import java.util.*;

public abstract class Node<T extends Region>
{
    private Node<T> lower = null, higher = null;
    public T value;

    public Node(T value)
    {
        this.value = value;
    }

    public Node(Iterable<T> values)
    {
        addAllOptimized(values);
    }

    public boolean matches(T term)
    {
        return value.equals(term);
    }

    public boolean contains(T value)
    {
        return getMatchingChild(value) != null;
    }

    public T getMatchingChild(T searchValue)
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

    public abstract void add(T value);

    public abstract void add(T value, boolean merge);

    public void add(Node<T> node)
    {
        add(node, false);
    }

    public void add(Node<T> node, boolean merge)
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

    public void addAllOptimized(Iterable<T> values)
    {
        addAllOptimized(values, false);
    }

    public void addAllOptimized(Iterable<T> values, boolean merge)
    {
        HashSet<T> valuesSet = new HashSet<>();
        values.forEach(valuesSet::add);
        valuesSet.addAll(getAllValues());

        ArrayList<T> valuesList = new ArrayList<>(valuesSet);

        Collections.sort(valuesList);

        addOptimizedRecursive(valuesList, merge);
    }

    private void addOptimizedRecursive(List<T> elements, boolean merge)
    {
        if (elements.isEmpty())
        {
            return;
        }

        int center = elements.size() / 2;

        List<T> lower = elements.subList(0, center);
        List<T> higher = elements.subList(center + 1, elements.size());

        add(elements.get(center), merge);
        addOptimizedRecursive(lower, merge);
        addOptimizedRecursive(higher, merge);
    }

    public List<T> getAllValuesSorted()
    {
        List<T> list = new ArrayList<>(getAllValues());
        Collections.sort(list);
        return list;
    }

    public Set<T> getAllValues()
    {
        Set<T> values = new HashSet<>();
        getAllValuesRecursive(values, this);
        return values;
    }

    private void getAllValuesRecursive(Set<T> values, Node<T> parent)
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
