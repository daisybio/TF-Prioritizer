package util.Regions;

/**
 * Stores the score of a peak as well as its corresponding region
 */
public class PeakRegion extends Region implements Comparable
{
    /**
     * The peak score
     */
    private final int score;
    /**
     * The line which was used in order to build this PeakRegion. Only set if the according constructor is used.
     */
    private String line;

    /**
     * Create a new PeakRegion based on its attributes.
     *
     * @param chromosome the chromosome the PeakRegion is located on
     * @param start      the start position of the region
     * @param end        the end position of the region
     * @param score      the peak score of this peak
     */
    public PeakRegion(String chromosome, int start, int end, int score)
    {
        super(chromosome, start, end);
        this.score = score;
    }

    /**
     * Create a new PeakRegion based on a line originating from a file
     *
     * @param line the source string
     */
    public PeakRegion(String line)
    {
        this(line.split("\t"));
        this.line = line;
    }

    /**
     * Create a new PeakRegion based on a splitted input line.
     *
     * @param split the split of the input line. Has to be at least 5 elements long.
     */
    private PeakRegion(String[] split)
    {
        super(split[0], Integer.parseInt(split[1]), Integer.parseInt(split[2]));

        if (split.length == 4)
        {
            this.score = Integer.parseInt(split[3]);
        } else
        {
            this.score = Integer.parseInt(split[4]);
        }
    }

    /**
     * @return get the score of this PeakRegion
     */
    public int getScore()
    {
        return score;
    }

    public String toString()
    {
        if (line != null)
        {
            return line;
        }
        return getChromosome() + "\t" + getStart() + "\t" + getEnd() + "\t" + getScore();
    }
}
