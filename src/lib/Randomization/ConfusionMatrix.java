package lib.Randomization;

public class ConfusionMatrix
{
    private int tn, tp, fp, fn;

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
}
