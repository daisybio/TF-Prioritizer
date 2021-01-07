package util;

public class BL_ranges_binary_tree implements Comparable
{
    public int number;
    public String chr;
    public int left_border;
    public int right_border;
    public String signal;

    public String toString()
    {
        String ret ="";

        ret+=number;
        ret+="\t";
        ret+=chr;
        ret+="\t";
        ret+=left_border;
        ret+="\t";
        ret+=right_border;
        ret+="\t";
        ret+=signal;

        return  ret;
    }

    @Override
    public int compareTo(Object o) {
        int compare = ((BL_ranges_binary_tree)o).left_border;

        return this.left_border-compare;
    }
}
