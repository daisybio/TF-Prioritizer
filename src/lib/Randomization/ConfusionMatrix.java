package lib.Randomization;

public class ConfusionMatrix
{
    private int tn;
    private int tp;
    private int fp;
    private int fn;
    private int tnRandomized;
    private int fnRandomized;

    public int getTnRandomized()
    {
        return tnRandomized;
    }

    public void setTnRandomized(int tnRandomized)
    {
        this.tnRandomized = tnRandomized;
    }

    public int getFnRandomized()
    {
        return fnRandomized;
    }

    public void setFnRandomized(int fnRandomized)
    {
        this.fnRandomized = fnRandomized;
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

    @Override public String toString()
    {
        return "ConfusionMatrix{" + "tn=" + tn + ", tp=" + tp + ", fp=" + fp + ", fn=" + fn + ", tnRandomized=" +
                tnRandomized + ", fnRandomized=" + fnRandomized + '}';
    }
}
