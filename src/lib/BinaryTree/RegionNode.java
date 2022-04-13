package lib.BinaryTree;

import lib.Region;
import util.BinarySearchTree.Node;

public class RegionNode extends Node<Region>
{
    public RegionNode(Region region)
    {
        super(region);
    }

    public RegionNode(Iterable<Region> regions)
    {
        super(regions);
    }

    @Override public boolean matches(Region term)
    {
        return !(Math.max(term.getStart(), term.getEnd()) < Math.min(value.getStart(), value.getEnd()) ||
                Math.min(term.getStart(), term.getEnd()) > Math.max(value.getStart(), value.getEnd()));
    }

    @Override public void add(Region value)
    {
        super.add(new RegionNode(value));
    }
}
