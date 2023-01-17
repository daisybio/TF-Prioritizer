package org.exbio.tfprio.steps.report;

import org.json.JSONObject;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

public class TfGroup {
    private final String symbol;
    private final Collection<TranscriptionFactor> transcriptionFactors;
    private final Map<String, Integer> hmRanks = new HashMap<>();

    public TfGroup(String symbol, Collection<TranscriptionFactor> transcriptionFactors) {
        this.symbol = symbol;
        this.transcriptionFactors = transcriptionFactors;
    }

    public void setRank(String hm, int rank) {
        hmRanks.put(hm, rank);
    }

    public Collection<TranscriptionFactor> getTranscriptionFactors() {
        return transcriptionFactors;
    }

    public String getSymbol() {
        return symbol;
    }

    public JSONObject toJSON() {
        return new JSONObject() {{
            put("symbol", symbol);
            put("hmRanks", hmRanks);
            put("transcriptionFactors", transcriptionFactors.stream().map(TranscriptionFactor::toJSON).toList());
        }};
    }
}
