package util.BinarySearchTree;

import lib.Region;

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
        super.add(new RegionNode(value), false);
    }

    @Override public void add(Region value, boolean merge)
    {
        super.add(new RegionNode(value), merge);
    }
}