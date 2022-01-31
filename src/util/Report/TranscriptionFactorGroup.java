package util.Report;

import java.util.ArrayList;
import java.util.Map;

public class TranscriptionFactorGroup
{
    private boolean hasValidation = true, hasDistribution = true, hasRegression = true;
    private final String name;
    private final ArrayList<TranscriptionFactor> transcriptionFactors;
    Map<String, Map<String, Number>> regressionCoefficients;
    boolean realGroup;

    TranscriptionFactorGroup(String name, ArrayList<TranscriptionFactor> transcriptionFactors,
                             Map<String, Map<String, Number>> regressionCoefficients)
    {
        this.name = name;
        this.transcriptionFactors = transcriptionFactors;
        this.regressionCoefficients = regressionCoefficients;
        this.realGroup = transcriptionFactors.size() > 1;
    }

    public boolean hasDistribution()
    {
        return hasDistribution;
    }

    public void setDistribution(boolean hasDistribution)
    {
        this.hasDistribution = hasDistribution;
    }

    public boolean hasValidation()
    {
        return hasValidation;
    }

    public void setValidation(boolean hasValidation)
    {
        this.hasValidation = hasValidation;
    }

    public boolean hasRegression()
    {
        return hasRegression;
    }

    public void setRegression(boolean hasRegression)
    {
        this.hasRegression = hasRegression;
    }

    public String getName()
    {
        return name;
    }

    public ArrayList<TranscriptionFactor> getTranscriptionFactors()
    {
        return transcriptionFactors;
    }

    public Map<String, Map<String, Number>> getRegressionCoefficients()
    {
        return regressionCoefficients;
    }
}