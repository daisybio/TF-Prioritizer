package lib.Blacklist;

import lib.Region;

public class BlacklistedRegion extends Region
{
    protected final String signal;
    protected Double score;

    public BlacklistedRegion(String chromosome, int start, int end, String signal)
    {
        super(chromosome, start, end);
        this.signal = signal;
    }

    public BlacklistedRegion(String[] split)
    {
        super(split[0], Integer.parseInt(split[1]), Integer.parseInt(split[2]));
        signal = split[3].replace(" ", "_").toUpperCase();
    }

    public String toString()
    {
        String output = "";

        output += getChromosome();
        output += "\t";
        output += getStart();
        output += "\t";
        output += getEnd();

        if (!signal.equals(""))
        {
            output += "\t";
            output += signal;
        }

        if (score != null)
        {
            output += "\t";
            output += score;
        }

        return output;
    }
}
