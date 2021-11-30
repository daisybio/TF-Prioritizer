package util;

public class Gene_Affinity_Value implements Comparable
{

    public String gene_name = "";
    public String gene_symbol = "NOT_AVAILABLE";
    public double affinity_value = 0.0;

    public String toString()
    {
        String res = "";
        res += gene_name;
        res += "\t";
        res += gene_symbol;
        res += "\t";
        res += affinity_value;

        return res;
    }


    @Override public int compareTo(Object o)
    {
        double compare = ((Gene_Affinity_Value) o).affinity_value;

        if (compare < this.affinity_value)
        {
            return -1;
        } else if (compare > this.affinity_value)
        {
            return 1;
        } else
        {
            return 0;
        }
    }
}
