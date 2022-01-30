package util.Report;

import java.util.ArrayList;
import java.util.Map;

public class TranscriptionFactor
{
    private final String geneID;
    private final String name;
    private final Map<String, Map<String, Number>> log2fc;
    private final Map<String, Number> tpm;
    private final Map<String, Number> normex;
    private final ArrayList<String> histoneModifications;

    public TranscriptionFactor(String geneID, String name, Map<String, Map<String, Number>> log2fc,
                               Map<String, Number> tpm, Map<String, Number> normex,
                               ArrayList<String> histoneModifications)
    {
        this.geneID = geneID;
        this.name = name;
        this.log2fc = log2fc;
        this.tpm = tpm;
        this.normex = normex;
        this.histoneModifications = histoneModifications;
    }

    public String getGeneID()
    {
        return geneID;
    }

    public String getName()
    {
        return name;
    }

    public Map<String, Map<String, Number>> getLog2fc()
    {
        return log2fc;
    }

    public Map<String, Number> getTpm()
    {
        return tpm;
    }

    public Map<String, Number> getNormex()
    {
        return normex;
    }

    public ArrayList<String> getHistoneModifications()
    {
        return histoneModifications;
    }
}
