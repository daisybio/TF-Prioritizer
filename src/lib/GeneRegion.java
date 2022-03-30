package lib;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class GeneRegion extends Region
{
    private final String id;

    public GeneRegion(String chromosome, int start, int end, String id)
    {
        super(chromosome, start, end);
        this.id = id.toUpperCase();
    }

    public String getId()
    {
        return id;
    }

    public void merge(GeneRegion other)
    {
        if (!getChromosome().equals(other.getChromosome()))
        {
            throw new IllegalArgumentException("Cannot merge regions on different chromosomes");
        }
        start = Math.min(start, other.getStart());
        end = Math.max(end, other.getEnd());
    }

    public boolean isIdentical(GeneRegion other, boolean idMatters)
    {
        return (!idMatters || other.getId().equals(getId())) && other.getChromosome().equals(getChromosome()) &&
                other.getStart() == getStart() && other.getEnd() == getEnd();
    }

    @Override public int compareTo(Object o)
    {
        GeneRegion region = (GeneRegion) o;

        int regionCompare = super.compareTo(region);
        if (regionCompare == 0)
        {
            return getId().compareTo(region.getId());
        } else
        {
            return regionCompare;
        }
    }

    public static List<GeneRegion> removeGeneRegionDuplicates(List<GeneRegion> input, boolean idMatters)
    {
        Collections.sort(input);
        List<GeneRegion> output = new ArrayList<>();
        for (int i = 0; i < input.size(); i++)
        {
            GeneRegion region = input.get(i);
            for (int j = i + 1; j < input.size(); j++)
            {
                GeneRegion other = input.get(j);
                if (!other.isIdentical(region, idMatters))
                {
                    i = j - 1;
                    break;
                }
            }
            output.add(region);
        }
        return output;
    }

    public Region getRegion()
    {
        return new Region(getChromosome(), getStart(), getEnd());
    }

    public String toString()
    {
        return getChromosome() + "\t" + getStart() + "\t" + getEnd() + "\t" + getId();
    }
}
