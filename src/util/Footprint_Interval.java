package util;

public class Footprint_Interval implements Comparable
{

    public String chromosome = "";
    public int start;
    public int end;
    public String name = "";
    public int score;
    public String strand = "";
    public double signalValue;
    public double pValue;
    public double qValue;
    public String line = "";

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
    }

    public void make_line()
    {
        StringBuilder sb = new StringBuilder();
        sb.append(chromosome);
        sb.append("\t");
        sb.append(start);
        sb.append("\t");
        sb.append(end);
        sb.append("\t");
        sb.append(name);
        sb.append("\t");
        sb.append(score);
        sb.append("\t");
        sb.append(strand);
        sb.append("\t");
        sb.append(signalValue);
        sb.append("\t");
        sb.append(pValue);
        sb.append("\t");
        sb.append(qValue);

        line = sb.toString();
    }

    @Override public int compareTo(Object o)
    {
        int compare = ((Footprint_Interval) o).start;

        return this.start - compare;
    }
}
