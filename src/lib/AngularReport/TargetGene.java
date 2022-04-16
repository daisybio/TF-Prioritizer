package lib.AngularReport;

import org.json.JSONObject;

public class TargetGene
{
    private final String geneID;
    private final String symbol;

    public TargetGene(String geneID, String symbol)
    {
        this.geneID = geneID;
        this.symbol = symbol;
    }

    public JSONObject toJSONObject()
    {
        return new JSONObject()
        {{
            accumulate("geneID", geneID);
            accumulate("symbol", symbol);
        }};
    }
}
