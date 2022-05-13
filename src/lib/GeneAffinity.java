package lib;

import tfprio.tfprio.TFPRIO;

public class GeneAffinity implements Comparable<GeneAffinity>
{
    private final String geneID, geneSymbol;
    private final Double value;
    private String geneClass;

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

    public String getGeneClass()
    {
        return geneClass;
    }

    public void setGeneClass(String geneClass)
    {
        this.geneClass = geneClass;
    }

    public String getGeneID()
    {
        return geneID;
    }

    public String getGeneSymbol()
    {
        return geneSymbol;
    }

    public Double getValue()
    {
        return value;
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
