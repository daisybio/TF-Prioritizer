package lib.DistributionAnalysis.Classes;

import java.util.ArrayList;
import java.util.Collections;

public class StatsCollection
{
    private final ArrayList<Stats> tfs = new ArrayList<>();
    private final Stats background;

    public StatsCollection(Stats background)
    {
        this.background = background;
    }

    public void sort()
    {
        Collections.sort(tfs);
    }

    public ArrayList<Stats> getTfs()
    {
        return tfs;
    }

    public Stats getBackground()
    {
        return background;
    }

    public void add(Stats stats)
    {
        tfs.add(stats);
    }
}
