package lib.DistributionAnalysis.Classes;

import java.util.ArrayList;
import java.util.List;

public class TargetGeneDcg implements Comparable<TargetGeneDcg>
{
    public final String targetGene;
    public final List<Double> affinityValues = new ArrayList<>();

    public TargetGeneDcg(String targetGene)
    {
        this.targetGene = targetGene;
    }

    public void addAffinityValue(double value)
    {
        affinityValues.add(value);
    }

    public Double getAverageAffinityValue()
    {
        double sum = 0.0;

        for (double affinityValue : affinityValues)
        {
            sum += affinityValue;
        }

        return sum / affinityValues.size();
    }

    public String getTargetGene()
    {
        return targetGene;
    }

    @Override public int compareTo(TargetGeneDcg other)
    {
        return other.getAverageAffinityValue().compareTo(getAverageAffinityValue());
    }

    @Override public String toString()
    {
        return targetGene + "\t" + getAverageAffinityValue();
    }
}
