package util;

public class Analysis_distribution_stats implements Comparable{
    public String label = "";
    public double sum_all_values = 0;
    public double number_target_genes = 0;
    public double mean = 0;
    public double median = 0;
    public double quantile_95 = 0;
    public double quantile_99 = 0;

    public Analysis_distribution_stats(){}


    @Override
    public int compareTo(Object o) {
        Analysis_distribution_stats other = (Analysis_distribution_stats) o;

        if(other.median > this.median)
        {
            return -1;
        }
        if(other.median < this.median)
        {
            return 1;
        }

        if(other.quantile_95 > this.quantile_95)
        {
            return -1;
        }
        if(other.quantile_95 < this.quantile_95)
        {
            return 1;
        }

        if(other.quantile_99 > this.quantile_99)
        {
            return -1;
        }
        else
        {
            return 1;
        }
    }

    public String to_html_String(Analysis_distribution_stats background)
    {
        StringBuilder sb = new StringBuilder();
        sb.append("\t\t<table style=\"width:80%\">\n");
        sb.append("\t\t\t<tr>\n");
        sb.append("\t\t\t\t<th>");
        sb.append("Distribution");
        sb.append("\t\t\t\t</th>\n");
        sb.append("\t\t\t\t<th>");
        sb.append("Mean");
        sb.append("\t\t\t\t</th>\n");
        sb.append("\t\t\t\t<th>");
        sb.append("Median");
        sb.append("\t\t\t\t</th>\n");
        sb.append("\t\t\t\t<th>");
        sb.append("95% Quantile");
        sb.append("\t\t\t\t</th>\n");
        sb.append("\t\t\t\t<th>");
        sb.append("99% Quantile");
        sb.append("\t\t\t\t</th>\n");
        sb.append("\t\t\t</tr>\n");
        sb.append("\t\t\t<tr>\n");
        sb.append("\t\t\t\t<th>");
        sb.append("BACKGROUND");
        sb.append("\t\t\t\t</th>\n");
        sb.append("\t\t\t\t<th>");
        sb.append(background.median);
        sb.append("\t\t\t\t</th>\n");
        sb.append("\t\t\t\t<th>");
        sb.append(background.mean);
        sb.append("\t\t\t\t</th>\n");
        sb.append("\t\t\t\t<th>");
        sb.append(background.quantile_95);
        sb.append("\t\t\t\t</th>\n");
        sb.append("\t\t\t\t<th>");
        sb.append(background.quantile_99);
        sb.append("\t\t\t\t</th>\n");
        sb.append("\t\t\t</tr>\n");
        sb.append("\t\t\t<tr>\n");
        sb.append("\t\t\t\t<th>");
        sb.append(this.label.toUpperCase());
        sb.append("\t\t\t\t</th>\n");
        sb.append("\t\t\t\t<th>");
        sb.append(this.median);
        sb.append("\t\t\t\t</th>\n");
        sb.append("\t\t\t\t<th>");
        sb.append(this.mean);
        sb.append("\t\t\t\t</th>\n");
        sb.append("\t\t\t\t<th>");
        sb.append(this.quantile_95);
        sb.append("\t\t\t\t</th>\n");
        sb.append("\t\t\t\t<th>");
        sb.append(this.quantile_99);
        sb.append("\t\t\t\t</th>\n");
        sb.append("\t\t\t</tr>\n");


        return  sb.toString();

    }
}
