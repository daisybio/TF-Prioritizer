package lib.ConfusionMatrixes;

import java.util.HashMap;
import java.util.Map;

public class ConfusionMatrix
{
    private int tn;
    private int tp;
    private int fp;
    private int fn;

    public static String getTabularHeader(String sep)
    {
        return String.join(sep, "TP", "TN", "FP", "FN", "Sensitivity", "Specificity", "Accuracy", "Precision", "F1");
    }

    public int getTn()
    {
        return tn;
    }

    public void setTn(int tn)
    {
        this.tn = tn;
    }

    public int getTp()
    {
        return tp;
    }

    public void setTp(int tp)
    {
        this.tp = tp;
    }

    public int getFp()
    {
        return fp;
    }

    public void setFp(int fp)
    {
        this.fp = fp;
    }

    public int getFn()
    {
        return fn;
    }

    public void setFn(int fn)
    {
        this.fn = fn;
    }

    private double getSpecificity()
    {
        return (double) getTn() / (getFp() + getTn());
    }

    private double getSensitivity()
    {
        return (double) getTp() / (getFn() + getTp());
    }

    private double getAccuracy()
    {
        return (double) (getTn() + getTp()) / (getTn() + getTp() + getFn() + getFp());
    }

    private double getPrecision()
    {
        return (double) getTp() / (getTp() + getFp());
    }

    private double getRecall()
    {
        return getSensitivity();
    }

    private double getF1()
    {
        return (double) 2 * (getPrecision() * getRecall()) / (getPrecision() + getRecall());
    }

    @Override public String toString()
    {
        return "ConfusionMatrix{" + "tn=" + tn + ", tp=" + tp + ", fp=" + fp + ", fn=" + fn + '}';
    }

    public Map<String, Double> toMap()
    {
        return new HashMap<>()
        {{
            put("tn", (double) tn);
            put("tp", (double) tp);
            put("fp", (double) fp);
            put("fn", (double) fn);
            put("sensitivity", getSensitivity());
            put("specificity", getSpecificity());
            put("accuracy", getAccuracy());
            put("precision", getPrecision());
            put("f1", getF1());
        }};
    }

    public String getTabular(String sep)
    {
        return String.join(sep, String.valueOf(getTp()), String.valueOf(getTn()), String.valueOf(getFp()),
                String.valueOf(getFn()), String.valueOf(getSensitivity()), String.valueOf(getSpecificity()),
                String.valueOf(getAccuracy()), String.valueOf(getPrecision()), String.valueOf(getF1()));
    }
}
