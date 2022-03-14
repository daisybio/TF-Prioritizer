package util;

import com2pose.COM2POSE;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

import static util.FileManagement.readFile;
import static util.FileManagement.writeFile;
import static util.ScriptExecution.executeAndWait;

public class MapSymbolAndEnsg
{
    private final Map<String, String> symbolEnsg = new HashMap<>();
    private final Map<String, String> ensgSymbol = new HashMap<>();
    private final Logger logger = new Logger(this.getClass().getSimpleName());

    private final File map = COM2POSE.configs.deSeq2.fileStructure.f_mapping.get();

    public MapSymbolAndEnsg()
    {
        try
        {
            if (!map.exists())
            {
                createMappingFile();
            }
            loadMappingFile();
        } catch (IOException | InterruptedException e)
        {
            e.printStackTrace();
        }
    }

    private void createMappingFile() throws IOException, InterruptedException
    {
        String script = readFile(COM2POSE.configs.scriptTemplates.f_mapping.get());

        script = script.replace("{INPUTFILE}", COM2POSE.configs.deSeq2.inputGeneID.get());
        script = script.replace("{DATASET_SPECIES}", COM2POSE.configs.deSeq2.biomartDatasetSpecies.get());
        script = script.replace("{SYMBOL_COLUMN}", COM2POSE.configs.deSeq2.biomartDatasetSymbolColumn.get());
        script =
                script.replace("{OUTPUTFILE}", COM2POSE.configs.deSeq2.fileStructure.f_mapping.get().getAbsolutePath());

        writeFile(COM2POSE.configs.deSeq2.fileStructure.f_mappingScript.get(), script);

        executeAndWait(COM2POSE.configs.deSeq2.fileStructure.f_mappingScript.get(), logger);
    }

    private void loadMappingFile() throws FileNotFoundException
    {
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
    }

    public String symbolToEnsg(String symbol) throws FileNotFoundException, NoSuchFieldException
    {
        if (symbolEnsg.containsKey(symbol.toUpperCase()))
        {
            return symbolEnsg.get(symbol.toUpperCase());
        } else
        {
            throw new NoSuchFieldException("Could not find symbol " + symbol + " in the map file.");
        }
    }

    public String ensgToSymbol(String ensg) throws FileNotFoundException, NoSuchFieldException
    {
        if (ensgSymbol.containsKey(ensg.toUpperCase()))
        {
            return ensgSymbol.get(ensg.toUpperCase());
        } else
        {
            throw new NoSuchFieldException("Could not find ensg " + ensg + " in the map file.");
        }
    }
}
