package lib.DistributionAnalysis.Classes;

import java.util.*;

public class TfTargetGeneDcg
{
    private final String tfSymbol;
    private final Map<String, TargetGeneDcg> targetGene_affinityValues = new HashMap<>();

    public TfTargetGeneDcg(String tfSymbol)
    {
        this.tfSymbol = tfSymbol;
    }

    public String getTfSymbol()
    {
        return tfSymbol;
    }

    public void addAffinityValue(String targetGene, double value)
    {
        if (!targetGene_affinityValues.containsKey(targetGene))
        {
            targetGene_affinityValues.put(targetGene, new TargetGeneDcg(targetGene));
        }
        targetGene_affinityValues.get(targetGene).addAffinityValue(value);
    }

    public boolean isEmpty()
    {
        return targetGene_affinityValues.isEmpty();
    }

    public List<TargetGeneDcg> getOrderedTargetGeneList()
    {
        List<TargetGeneDcg> output = new ArrayList<>();

        for (String targetGene : targetGene_affinityValues.keySet())
        {
            output.add(targetGene_affinityValues.get(targetGene));
        }

        Collections.sort(output);

        return output;
    }
}
