package util.Regions;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Stores the id of a gene as well as its corresponding region
 */
public class GeneRegion extends Region
{
    /**
     * The gene id of this GeneRegion
     */
    private final String id;

    /**
     * Creat a new GeneRegion based on its attributes
     *
     * @param chromosome the chromosome this region is located on
     * @param start      the start position of this region
     * @param end        the end position of this region
     * @param id         the gene id of this gene
     */
    public GeneRegion(String chromosome, int start, int end, String id)
    {
        super(chromosome, start, end);
        this.id = id.toUpperCase();
    }

    /**
     * @return the gene id of this gene
     */
    public String getId()
    {
        return id;
    }

    /**
     * Check if this GeneRegion is identical with a given other GeneRegion
     *
     * @param other     the GeneRegion to compare this GeneRegion to
     * @param idMatters indicates if id should be considered
     * @return true if all the considered attributes are equal, otherwise false
     */
    public boolean isIdentical(GeneRegion other, boolean idMatters)
    {
        return (!idMatters || other.getId().equals(getId())) && other.getChromosome().equals(getChromosome()) &&
                other.getStart() == getStart() && other.getEnd() == getEnd();
    }

    /**
     * Compare this GeneRegion to a given other GeneRegion
     * <p>
     * The comparison is made in the following order:
     * <ol>
     *     <li>Chromosome</li>
     *     <li>Region start position</li>
     *     <li>Gene id</li>
     * </ol>
     *
     * @param o the object to be compared.
     * @return the comparison result
     */
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

    /**
     * Remove all duplicate geneRegions from a given list of GeneRegions.
     *
     * @param input     the list to remove duplicates from
     * @param idMatters indicates if id should be considered when checking for duplicates
     * @return the filtered list of geneRegions
     */
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

    /**
     * Get a {@link Region} object based on the attributes of this object
     *
     * @return the newly created Region object
     */
    public Region getRegion()
    {
        return new Region(getChromosome(), getStart(), getEnd());
    }

    public String toString()
    {
        return getChromosome() + "\t" + getStart() + "\t" + getEnd() + "\t" + getId();
    }

    public static List<MultiGeneRegion> mergeSameRegions(List<GeneRegion> input)
    {
        List<MultiGeneRegion> combined = new ArrayList<>();

        for (int i = 0; i < input.size(); i++)
        {
            GeneRegion first = input.get(i);
            MultiGeneRegion multiGeneRegion =
                    new MultiGeneRegion(first.getChromosome(), first.getStart(), first.getEnd(), first.getId());

            for (int j = i + 1; j < input.size(); j++)
            {
                GeneRegion current = input.get(j);

                if (first.isIdentical(current, false))
                {
                    multiGeneRegion.addId(current.getId());
                    i++;
                } else
                {
                    break;
                }
            }

            combined.add(multiGeneRegion);
        }
        return combined;
    }
}
