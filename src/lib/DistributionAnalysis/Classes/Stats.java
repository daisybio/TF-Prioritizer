package lib.DistributionAnalysis.Classes;

public class Stats implements Comparable<Stats>
{
    private final String label;
    private final double sum_all_values;
    private final double number_target_genes;
    private final double mean;
    private final Double median;
    private final Double quantile_95;
    private final Double quantile_99;

    public Stats(String line)
    {
        String[] split = line.split("\t");
        label = split[1];
        sum_all_values = Double.parseDouble(split[2]);
        number_target_genes = Double.parseDouble(split[3]);
        mean = Double.parseDouble(split[4]);
        median = Double.parseDouble(split[5]);
        quantile_95 = Double.parseDouble(split[6]);
        quantile_99 = Double.parseDouble(split[7]);
    }

    public Stats(String label, double sum_all_values, double number_target_genes, double mean, double median,
                 double quantile_95, double quantile_99)
    {
        this.label = label;
        this.sum_all_values = sum_all_values;
        this.number_target_genes = number_target_genes;
        this.mean = mean;
        this.median = median;
        this.quantile_95 = quantile_95;
        this.quantile_99 = quantile_99;
    }

    public String getLabel()
    {
        return label;
    }

    public double getSum_all_values()
    {
        return sum_all_values;
    }

    public double getNumber_target_genes()
    {
        return number_target_genes;
    }

    public double getMean()
    {
        return mean;
    }

    public Double getMedian()
    {
        return median;
    }

    public Double getQuantile_95()
    {
        return quantile_95;
    }

    public Double getQuantile_99()
    {
        return quantile_99;
    }

    @Override public int compareTo(Stats other)
    {
        if (other.getMedian().compareTo(getMedian()) != 0)
        {
            return other.getMedian().compareTo(getMedian());
        }

        if (other.getQuantile_95().compareTo(getQuantile_95()) != 0)
        {
            return other.getQuantile_95().compareTo(getQuantile_95());
        }

        return other.getQuantile_99().compareTo(getQuantile_99());
    }
}
