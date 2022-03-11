package util.BinarySearchTree;

public abstract class Node<T extends Comparable<T>>
{
    private Node<T> lower = null, higher = null;
    public final T value;

    public Node(T value)
    {
        this.value = value;
    }

    public boolean matches(T term)
    {
        return value.equals(term);
    }

    public boolean contains(T searchValue)
    {
        if (matches(searchValue))
        {
            return true;
        } else
        {
            int compareValue = value.compareTo(searchValue);
            if (compareValue < 0)
            {
                if (higher != null)
                {
                    return higher.contains(searchValue);
                } else
                {
                    return false;
                }
            } else if (compareValue > 0)
            {
                if (lower != null)
                {
                    return lower.contains(searchValue);
                } else
                {
                    return false;
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

    public String toString()
    {
        return toString(0);
    }

    public String toString(int level)
    {
        String output = (higher != null) ? higher.toString(level + 1) : "\t".repeat(level + 1) + "NULL\n";
        output += "\t".repeat(level) + this.value + "\n";
        output += (lower != null) ? lower.toString(level + 1) : "\t".repeat(level + 1) + "NULL\n";
        return output;
    }
}
