package lib.MixOptions;

public class MIX_Interval_Object
{
    String chr;
    public int start;
    public int end;
    String description;
    int val1;
    String val2;
    double val3;
    double val4;
    double val5;

    public MIX_Interval_Object(String chr, int start, int end, String descriptions, int val1, String val2, double val3,
                               double val4, double val5)
    {
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.description = descriptions;
        this.val1 = val1;
        this.val2 = val2;
        this.val3 = val3;
        this.val4 = val4;
        this.val5 = val5;
    }

    public MIX_Interval_Object()
    {
    }

    ;
}
