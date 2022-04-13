package lib;

public class Peak extends Region implements Comparable
{
    private final int score;
    private String line;

    public Peak(String chromosome, int start, int end, int score)
    {
        super(chromosome, start, end);
        this.score = score;
    }

    public Peak(String line)
    {
        this(line.split("\t"));
        this.line = line;
    }

    private Peak(String[] split)
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
