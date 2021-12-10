package util;

public class DESeq2_gene_object implements Comparable
{
    public String ensg_name ="";
    public double log2fc = 0;


    @Override public int compareTo(Object o)
    {
        double compare = ((DESeq2_gene_object) o).log2fc;

        if (compare < this.log2fc)
        {
            return -1;
        } else if (compare > this.log2fc)
        {
            return 1;
        } else
        {
            return 0;
        }
    }
}
