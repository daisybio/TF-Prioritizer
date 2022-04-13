package util.FileFilters;

public class Filters
{
    public static DirectoryFilter directoryFilter = new DirectoryFilter();
    public static FileFilter fileFilter = new FileFilter();

    public static java.io.FileFilter getSuffixFilter(String suffix)
    {
        return new SuffixFilter(suffix);
    }
}
