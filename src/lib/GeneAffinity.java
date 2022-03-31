package lib;

import tfprio.TFPRIO;

public class GeneAffinity implements Comparable<GeneAffinity>
{
    private final String geneID, geneSymbol;
    private final Double value;

    public GeneAffinity(String geneID, double value)
    {
        this.geneID = geneID;
        this.value = value;

        if (TFPRIO.mapSymbolAndEnsg.hasEnsg(geneID))
        {
            String tmp = null;
            try
            {
                tmp = TFPRIO.mapSymbolAndEnsg.ensgToSymbol(geneID);
            } catch (NoSuchFieldException ignore)
            {
            } finally
            {
                geneSymbol = tmp == null ? "NO_SYMBOL" : tmp;
            }
        } else
        {
            geneSymbol = "NO_SYMBOL";
        }
    }

    @Override public String toString()
    {
        return geneID + "\t" + geneSymbol + "\t" + value;
    }

    @Override public int compareTo(GeneAffinity other)
    {
        return other.value.compareTo(value);
    }
}
