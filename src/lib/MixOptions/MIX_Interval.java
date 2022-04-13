package lib.MixOptions;

import java.util.ArrayList;

public class MIX_Interval implements Comparable<MIX_Interval>
{
    public int start;
    public int end;

    public ArrayList<MIX_Interval_Object> merged_intervals = new ArrayList<>();

    public MIX_Interval_Object means;

    public MIX_Interval(int start, int end)
    {
        this.start = start;
        this.end = end;
    }

    public void calculate_mean(String level)
    {
        means = new MIX_Interval_Object();

        means.start = start;
        means.end = end;
        means.chr = merged_intervals.get(0).chr;

        String description_pre = merged_intervals.get(0).description;
        if (level.equals("SAMPLE_LEVEL"))
        {
            String[] split = description_pre.split("_");
            if (split.length > 4)
            {
                description_pre = split[0] + "_" + split[1] + "_" + split[2] + "_SMIX_" + split[4] + "_";
            } else
            {
                description_pre = description_pre;
            }
        }
        if (level.equals("HM_LEVEL"))
        {
            String[] split = description_pre.split("_");
            if (split.length > 4)
            {
                description_pre = split[0] + "_" + split[1] + "_HMMIX_" + split[3] + "_" + split[4] + "_";
            } else
            {
                description_pre = description_pre;
            }
        }
        means.description = description_pre;
        means.val2 = merged_intervals.get(0).val2;

        int val1_mean = 0;
        double val3_mean = 0;
        double val4_mean = 0;
        double val5_mean = 0;

        for (MIX_Interval_Object merged_interval : merged_intervals)
        {
            val1_mean += merged_interval.val1;
            val3_mean += merged_interval.val3;
            val4_mean += merged_interval.val4;
            val5_mean += merged_interval.val5;
        }

        val1_mean /= merged_intervals.size();
        val3_mean /= merged_intervals.size();
        val4_mean /= merged_intervals.size();
        val5_mean /= merged_intervals.size();

        means.val1 = val1_mean;
        means.val3 = val3_mean;
        means.val4 = val4_mean;
        means.val5 = val5_mean;
    }

    public String meanToString(int peak)
    {
        return means.chr + "\t" + means.start + "\t" + means.end + "\t" + means.description + peak + "\t" + means.val1 +
                "\t" + means.val2 + "\t" + means.val3 + "\t" + means.val4 + "\t" + means.val5;
    }

    @Override public int compareTo(MIX_Interval compare)
    {
        return this.start - compare.start;
    }
}
