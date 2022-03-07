package util.FileFilters;

import java.io.File;

public class FileFilter implements java.io.FileFilter
{
    @Override public boolean accept(File file)
    {
        return file.isFile();
    }
}
