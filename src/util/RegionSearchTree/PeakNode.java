package util.RegionSearchTree;

import util.Regions.PeakRegion;
import util.Regions.Region;

import java.util.Set;

/**
 * A binary tree node which represents a certain GeneRegion
 */
public class PeakNode extends RegionNode
{
    /**
     * Create a new node based on a single peakRegion
     *
     * @param data the peakRegion
     */
    public PeakNode(Region data)
    {
        super(data);
    }

    /**
     * Create a new node structure based on multiple PeakRegions
     *
     * @param peaks the peakRegions to add
     */
    public PeakNode(Iterable<Region> peaks)
    {
        super(peaks);
    }

    /**
     * @return the average peak score of this node and all its child nodes
     */
    public double getAveragePeakScore()
    {
        Set<Region> peaks = getAllValues();

        int sum = 0;

        for (Region region : peaks)
        {
            sum += ((PeakRegion) region).getScore();
        }
        return (double) sum / peaks.size();
    }
}
