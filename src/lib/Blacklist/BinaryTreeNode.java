package lib.Blacklist;

import lib.Region;
import util.BinarySearchTree.Node;

public class BinaryTreeNode extends Node<Region>
{
    public BinaryTreeNode(Region data)
    {
        super(data);
    }


    @Override public boolean matches(Region term)
    {
        return !(Math.max(term.getStart(), term.getEnd()) < Math.min(value.getStart(), value.getEnd()) ||
                Math.min(term.getStart(), term.getEnd()) > Math.max(value.getStart(), value.getEnd()));
    }

    @Override public void add(Region value)
    {
        super.add(new BinaryTreeNode(value));
    }
}
