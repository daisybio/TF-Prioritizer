package util;

import tfprio.TFPRIO;
import util.Configs.Config;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
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

    private final Config<File> f_map = TFPRIO.configs.deSeq2.fileStructure.f_mapping;

    public MapSymbolAndEnsg()
    {
        try
        {
            if (!f_map.get().exists())
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
        String script = readFile(TFPRIO.configs.scriptTemplates.f_mapping.get());
        Config<File> f_script = TFPRIO.configs.deSeq2.fileStructure.f_mappingScript;

        script = script.replace("{INPUTFILE}", TFPRIO.configs.deSeq2.inputGeneID.get().getAbsolutePath());
        script = script.replace("{DATASET_SPECIES}", TFPRIO.configs.deSeq2.biomartDatasetSpecies.get());
        script = script.replace("{SYMBOL_COLUMN}", TFPRIO.configs.deSeq2.biomartDatasetSymbolColumn.get());
        script = script.replace("{OUTPUTFILE}", f_map.get().getAbsolutePath());

        writeFile(f_script.get(), script);

        TFPRIO.createdFileStructure.addAll(Arrays.asList(f_script, f_map));

        executeAndWait(TFPRIO.configs.deSeq2.fileStructure.f_mappingScript.get(), logger);
    }

    private void loadMappingFile() throws FileNotFoundException
    {
        try (Scanner scanner = new Scanner(f_map.get()))
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

    public String symbolToEnsg(String symbol) throws NoSuchFieldException
    {
        if (symbolEnsg.containsKey(symbol.toUpperCase()))
        {
            return symbolEnsg.get(symbol.toUpperCase());
        } else
        {
            throw new NoSuchFieldException("Could not find symbol " + symbol + " in the map file.");
        }
    }

    public String ensgToSymbol(String ensg) throws NoSuchFieldException
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
