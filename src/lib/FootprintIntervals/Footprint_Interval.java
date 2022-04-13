package lib.FootprintIntervals;

public class Footprint_Interval implements Comparable<Footprint_Interval>
{

    final String chromosome;
    final int start;
    final int end;
    final String name;
    final int score;
    final String strand;
    final double signalValue;
    final double pValue;
    final double qValue;
    private final String line;

    public Footprint_Interval(String chromosome, int start, int end, String name, int score, String strand,
                              double signalValue, double pValue, double qValue, String line)
    {
        this.chromosome = chromosome;
        this.start = start;
        this.end = end;
        this.name = name;
        this.score = score;
        this.strand = strand;
        this.signalValue = signalValue;
        this.pValue = pValue;
        this.qValue = qValue;
        this.line = line;
    }

    public Footprint_Interval(String chromosome, int start, int end, String name, int score, String strand,
                              double signalValue, double pValue, double qValue)
    {
        this.chromosome = chromosome;
        this.start = start;
        this.end = end;
        this.name = name;
        this.score = score;
        this.strand = strand;
        this.signalValue = signalValue;
        this.pValue = pValue;
        this.qValue = qValue;
        this.line = chromosome + "\t" + start + "\t" + end + "\t" + name + "\t" + score + "\t" + strand + "\t" +
                signalValue + "\t" + pValue + "\t" + qValue;
    }

    public String toString()
    {
        return line;
    }

    @Override public int compareTo(Footprint_Interval o)
    {
        return this.start - o.start;
    }
}
