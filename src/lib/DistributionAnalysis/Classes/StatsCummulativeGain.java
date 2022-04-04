package lib.DistributionAnalysis.Classes;

import java.util.Map;

public class StatsCummulativeGain implements Comparable<StatsCummulativeGain>
{
    private final String label;
    private final double score;
    private final Map<String, Stats> group_object;
    private final Map<String, Stats> group_background;

    public StatsCummulativeGain(String label, double score, Map<String, Stats> group_object,
                                Map<String, Stats> group_background)
    {
        this.label = label;
        this.score = score;
        this.group_object = group_object;
        this.group_background = group_background;
    }

    @Override public String toString()
    {
        return label + '\t' + score;
    }

    public String getLabel()
    {
        return label;
    }

    public Double getScore()
    {
        return score;
    }

    public Map<String, Stats> getGroup_object()
    {
        return group_object;
    }

    public Map<String, Stats> getGroup_background()
    {
        return group_background;
    }

    @Override public int compareTo(StatsCummulativeGain other)
    {
        return other.getScore().compareTo(getScore());
    }
}
