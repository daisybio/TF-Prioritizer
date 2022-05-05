package util.Regions;

import util.Comparators.ChromosomeComparator;

/**
 * Represents a certain region on a certain chromosome.
 */
public class Region implements Comparable
{
    /**
     * The chromosome this region is located on.
     */
    private final String chromosome;

    /**
     * The start position of this region
     */
    protected int start;
    /**
     * The end position of this region
     */
    protected int end;

    /**
     * Creates a new Region based on the given input data
     *
     * @param chromosome the chromosome this region is based on
     * @param start      the start position of this region
     * @param end        the end position of this region
     */
    public Region(String chromosome, int start, int end)
    {
        this.chromosome = chromosome.replace("chr", "");

        // Make sure start position is smaller than end position
        this.start = Math.min(start, end);
        this.end = Math.max(start, end);
    }

    /**
     * Create a new region based on an array of strings.
     * <p>
     * Is used if a new region is created during file reading processes.
     *
     * @param split the base string array. Has to be at least 3 entries long. The order of the entries has to be like
     *              this: <ol>
     *              <li>Chromosome</li>
     *              <li>Start position</li>
     *              <li>End position</li>
     *              </ol>
     */
    public Region(String[] split)
    {
        this(split[0], Integer.parseInt(split[1]), Integer.parseInt(split[2]));
    }

    /**
     * @return the chromosome this region is located on
     */
    public String getChromosome()
    {
        return chromosome;
    }

    /**
     * @return the start position of this region
     */
    public int getStart()
    {
        return start;
    }

    /**
     * @return the end position of this region
     */
    public int getEnd()
    {
        return end;
    }

    /**
     * Check if this region overlaps a given other region.
     *
     * @param other the region to check for overlaps
     * @return true if the regions overlap, otherwise false
     */
    public boolean overlaps(Region other)
    {
        if (getChromosome().equals(other.getChromosome()))
        {
            return Math.max(other.getStart(), other.getEnd()) >= Math.min(getStart(), getEnd()) &&
                    Math.min(other.getStart(), other.getEnd()) <= Math.max(getStart(), getEnd());
        }
        return false;
    }

    /**
     * Compare this region to a given other region.
     * <p>
     * First the chromosomes are considered, then the start position on the chromosome.
     *
     * @param o the object to be compared.
     * @return the comparison result
     */
    @Override public int compareTo(Object o)
    {
        Region region = (Region) o;

        if (!this.getChromosome().equals(region.getChromosome()))
        {
            return new ChromosomeComparator().compare(this.getChromosome(), region.getChromosome());
        }
        return this.getStart() - region.getStart();
    }

    /**
     * Check if a region has the same chromosome, start position and end position as a given other region
     *
     * @param other the region to check for identity
     * @return true if all the attributes are equal, otherwise false
     */
    public boolean isIdentical(Region other)
    {
        return other.getChromosome().equals(getChromosome()) && other.getStart() == getStart() &&
                other.getEnd() == getEnd();
    }

    /**
     * Extend the borders of this regions so that its original span as well as the span of the added region and
     * everything between them is covered.
     *
     * @param other the region to merge into this region
     */
    public void merge(Region other)
    {
        if (!getChromosome().equals(other.getChromosome()))
        {
            throw new IllegalArgumentException("Trying to merge regions of different chromosomes");
        }

        this.start = Math.min(getStart(), other.getStart());
        this.end = Math.max(getEnd(), other.getEnd());
    }
}
