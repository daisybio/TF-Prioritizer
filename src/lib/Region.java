package lib;

import util.Comparators.ChromosomeComparator;

public class Region implements Comparable
{
    private final String chromosome;
    protected int start, end;

    public Region(String chromosome, int start, int end)
    {
        this.chromosome = chromosome.replace("chr", "");
        this.start = Math.min(start, end);
        this.end = Math.max(start, end);
    }

    public Region(String[] split)
    {
        this(split[0], Integer.parseInt(split[1]), Integer.parseInt(split[2]));
    }

    public String getChromosome()
    {
        return chromosome;
    }

    public int getStart()
    {
        return start;
    }

    public int getEnd()
    {
        return end;
    }

    public boolean overlaps(Region other)
    {
        if (getChromosome().equals(other.getChromosome()))
        {
            return Math.max(other.getStart(), other.getEnd()) >= Math.min(getStart(), getEnd()) &&
                    Math.min(other.getStart(), other.getEnd()) <= Math.max(getStart(), getEnd());
        }
        return false;
    }

    @Override public int compareTo(Object o)
    {
        Region region = (Region) o;

        if (!this.getChromosome().equals(region.getChromosome()))
        {
            return new ChromosomeComparator().compare(this.getChromosome(), region.getChromosome());
        }
        return this.getStart() - region.getStart();
    }

    public boolean isIdentical(Region other)
    {
        return other.getChromosome().equals(getChromosome()) && other.getStart() == getStart() &&
                other.getEnd() == getEnd();
    }
}
