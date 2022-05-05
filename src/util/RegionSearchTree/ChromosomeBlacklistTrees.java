package util.RegionSearchTree;

import lib.Blacklist.BlacklistedRegion;
import lib.Region;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

public class ChromosomeBlacklistTrees extends ChromosomeRegionTrees
{
    public ChromosomeBlacklistTrees(File source)
    {
        addFile(source);
    }

    private void addFile(File source)
    {
        Set<Region> regions = new HashSet<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(source)))
        {
            String inputLine;
            while ((inputLine = reader.readLine()) != null)
            {
                regions.add(new BlacklistedRegion(inputLine));
            }
        } catch (IOException e)
        {
            e.printStackTrace();
        }
        addAllOptimized(regions);
    }
}
