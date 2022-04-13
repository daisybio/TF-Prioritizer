package util.BinarySearchTree;

import java.util.*;

public abstract class Node<T extends Comparable>
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

    public void add(Node<T> node)
    {
        if (this.value == null)
        {
            this.value = node.value;
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
        HashSet<T> valuesSet = new HashSet<>();
        values.forEach(valuesSet::add);
        valuesSet.addAll(getAllValues());

        ArrayList<T> valuesList = new ArrayList<>(valuesSet);

        Collections.sort(valuesList);

        addOptimizedRecursive(valuesList);
    }

    private void addOptimizedRecursive(List<T> elements)
    {
        if (elements.isEmpty())
        {
            return;
        }

        int center = elements.size() / 2;

        List<T> lower = elements.subList(0, center);
        List<T> higher = elements.subList(center + 1, elements.size());

        add(elements.get(center));
        addOptimizedRecursive(lower);
        addOptimizedRecursive(higher);
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
