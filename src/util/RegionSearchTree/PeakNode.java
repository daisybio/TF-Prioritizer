package util.RegionSearchTree;

import lib.Peak;
import lib.Region;

import java.util.Set;

public class PeakNode extends RegionNode
{
    public PeakNode(Region data)
    {
        super(data);
    }

    public PeakNode(Iterable<Region> peaks)
    {
        super(peaks);
    }

    public double getAveragePeakScore()
    {
        Set<Region> peaks = getAllValues();

        int sum = 0;

        for (Region region : peaks)
        {
            sum += ((Peak) region).getScore();
        }
        return (double) sum / peaks.size();
    }
}
