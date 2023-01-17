package org.exbio.tfprio.steps.report;

import org.json.JSONObject;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class TranscriptionFactor {
    private final List<String> ensgs;
    private final String symbol;
    private final Map<String, Double> pairingLog2fc = new HashMap<>();
    private final Map<String, Double> groupMeanExpression = new HashMap<>();
    private final Map<String, Double> groupTpm = new HashMap<>();

    public TranscriptionFactor(String symbol, List<String> ensgs) {
        this.symbol = symbol;
        this.ensgs = ensgs;
    }

    public void setLog2fc(String pairing, double log2fc) {
        this.pairingLog2fc.put(pairing, log2fc);
    }

    public void setMeanExpression(String group, double meanExpression) {
        this.groupMeanExpression.put(group, meanExpression);
    }

    public void setTpm(String group, double tpm) {
        this.groupTpm.put(group, tpm);
    }

    public String getSymbol() {
        return symbol;
    }

    public List<String> getEnsgs() {
        return ensgs;
    }

    public JSONObject toJSON() {
        return new JSONObject() {{
            put("symbol", symbol);
            put("ensgs", ensgs);
            put("pairingLog2fc", pairingLog2fc);
            put("groupMeanExpression", groupMeanExpression);
            put("groupTpm", groupTpm);
        }};
    }
}
