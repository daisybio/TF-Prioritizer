package util.Comparators;

public class ChromosomeComparator extends IntegerStringComparator
{
    @Override public int compare(String a, String b)
    {
        if (a.equals("MT") ^ b.equals("MT")) // exclusive or
        {
            return (a.equals("MT") ? 1 : -1);
        } else // if none equals MT or both equal MT
        {
            return super.compare(a, b);
        }
    }
}
