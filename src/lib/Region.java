package lib;

public class Region implements Comparable
{
    private final String chromosome;
    private final int start, end;

    public Region(String chromosome, int start, int end)
    {
        this.chromosome = chromosome.replace("chr", "");
        this.start = start;
        this.end = end;
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

    @Override public int compareTo(Object o)
    {
        Region region = (Region) o;

        if (!this.getChromosome().equals(region.getChromosome()))
        {
            throw new IllegalArgumentException("Trying to compare regions on different chromosomes");
        }
        return this.getStart() - region.getStart();
    }
}
