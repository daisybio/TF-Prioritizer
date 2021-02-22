package util;

import java.util.HashMap;

public class Analysis_distribution_stats_cumulative_gain implements Comparable {

    public String label;
    public int rank;
    public HashMap<String,Analysis_distribution_stats> group_object;
    public HashMap<String,Analysis_distribution_stats> group_background;

    @Override
    public int compareTo(Object o) {
        Analysis_distribution_stats_cumulative_gain other = (Analysis_distribution_stats_cumulative_gain) o;

        if(other.rank > this.rank)
        {
            return -1;
        }
        else
        {
            return 1;
        }
    }
}
