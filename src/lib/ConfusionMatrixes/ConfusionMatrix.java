package lib.ConfusionMatrixes;

import java.util.HashMap;
import java.util.Map;

public class ConfusionMatrix
{
    private int tn;
    private int tp;
    private int fp;
    private int fn;

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

    @Override public String toString()
    {
        return "ConfusionMatrix{" + "tn=" + tn + ", tp=" + tp + ", fp=" + fp + ", fn=" + fn + '}';
    }

    public Map<String, Integer> toMap()
    {
        return new HashMap<>()
        {{
            put("tn", tn);
            put("tp", tp);
            put("fp", fp);
            put("fn", fn);
        }};
    }
}
