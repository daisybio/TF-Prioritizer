package util.Comparators;

import java.io.File;
import java.util.Comparator;

public class FileComparator implements Comparator<File>
{
    @Override public int compare(File first, File second)
    {
        return first.getAbsolutePath().compareTo(second.getAbsolutePath());
    }
}
