package util;

public class BL_ranges_binary_tree implements Comparable {
    public int number;
    public String chr;
    public int left_border;
    public int right_border;
    public String signal = "";

    public double peak_score = -1.0;

    public String toString() {
        String ret = "";

        ret += number;
        ret += "\t";
        ret += chr;
        ret += "\t";
        ret += left_border;
        ret += "\t";
        ret += right_border;
        if (!signal.equals("")) {
            ret += "\t";
            ret += signal;
        }
        if (peak_score > -1) {
            ret += "\t";
            ret += peak_score;
        }

        return ret;
    }

    @Override public int compareTo(Object o) {
        int compare = ((BL_ranges_binary_tree) o).left_border;

        return this.left_border - compare;
    }
}
