package util;

import com2pose.COM2POSE;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

public class MapSymbolAndEnsg
{
    private static final Map<String, String> symbolEnsg = new HashMap<>();
    private static final Map<String, String> ensgSymbol = new HashMap<>();
    private static boolean initialized = false;

    private static void initialize() throws FileNotFoundException
    {
        File map = COM2POSE.configs.deSeq2.fileStructure.f_mapping.get();

        try (Scanner scanner = new Scanner(map))
        {
            while (scanner.hasNextLine())
            {
                String[] line = scanner.nextLine().split("\t");
                if (line.length > 1)
                {
                    symbolEnsg.put(line[1].toUpperCase(), line[0].toUpperCase());
                    ensgSymbol.put(line[0].toUpperCase(), line[1].toUpperCase());
                }
            }
        }

        initialized = true;
    }

    public static String symbolToEnsg(String symbol) throws FileNotFoundException, NoSuchFieldException
    {
        if (!initialized)
        {
            initialize();
        }
        if (symbolEnsg.containsKey(symbol.toUpperCase()))
        {
            return symbolEnsg.get(symbol.toUpperCase());
        } else
        {
            throw new NoSuchFieldException("Could not find symbol " + symbol + " in the map file.");
        }
    }

    public static String ensgToSymbol(String ensg) throws FileNotFoundException, NoSuchFieldException
    {
        if (!initialized)
        {
            initialize();
        }
        if (ensgSymbol.containsKey(ensg.toUpperCase()))
        {
            return ensgSymbol.get(ensg.toUpperCase());
        } else
        {
            throw new NoSuchFieldException("Could not find ensg " + ensg + " in the map file.");
        }
    }
}
