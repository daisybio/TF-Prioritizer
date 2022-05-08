package util;

import tfprio.tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

import static util.FileManagement.readFile;
import static util.FileManagement.writeFile;
import static util.ScriptExecution.executeAndWait;

public class MapSymbolAndEnsg
{
    private final Map<String, Set<String>> symbolEnsg = new HashMap<>();
    private final Map<String, String> ensgSymbol = new HashMap<>();
    private final Map<String, String> ensgDescription = new HashMap<>();
    private final Map<String, String> symbolDescription = new HashMap<>();
    private final Logger logger = new Logger(this.getClass().getSimpleName());
    private final List<String> ensgList = new ArrayList<>();

    private final AbstractConfig<File> f_map = TFPRIO.configs.deSeq2.fileStructure.f_mapping;
    private final AbstractConfig<File> f_scriptTemplate = TFPRIO.configs.scriptTemplates.f_mapping;
    private final AbstractConfig<File> f_script = TFPRIO.configs.deSeq2.fileStructure.f_mappingScript;
    private final AbstractConfig<File> f_geneIDs = TFPRIO.configs.deSeq2.inputGeneID;
    private final AbstractConfig<String> datasetSpecies = TFPRIO.configs.deSeq2.biomartDatasetSpecies;
    private final AbstractConfig<String> datasetSymbolColumn = TFPRIO.configs.deSeq2.biomartDatasetSymbolColumn;

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
        String script = readFile(f_scriptTemplate.get());

        script = script.replace("{INPUTFILE}", f_geneIDs.get().getAbsolutePath());
        script = script.replace("{DATASET_SPECIES}", datasetSpecies.get());
        script = script.replace("{SYMBOL_COLUMN}", datasetSymbolColumn.get());
        script = script.replace("{OUTPUTFILE}", f_map.get().getAbsolutePath());

        writeFile(f_script.get(), script);

        TFPRIO.createdFileStructure.addAll(Arrays.asList(f_script, f_map));

        executeAndWait(f_script.get(), logger);
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
                    String ensg = line[0].toUpperCase();
                    String symbol = line[1].toUpperCase();
                    symbolEnsg.put(symbol, new HashSet<>(List.of(ensg)));
                    ensgSymbol.put(ensg, symbol);
                    ensgList.add(ensg);
                    String description = line.length > 2 ? line[2] : "";
                    ensgDescription.put(ensg, description);
                    symbolDescription.put(symbol, description);

                    if (symbol.contains("-"))
                    {
                        String correctedSymbol = symbol.replace("-", ".");

                        symbolEnsg.put(correctedSymbol, new HashSet<>(List.of(ensg)));
                        symbolDescription.put(correctedSymbol, description);
                    }
                }
            }
        }
    }

    public Set<String> symbolToEnsg(String symbol) throws NoSuchFieldException
    {
        Map<String, String> groupSeparators = new HashMap<>()
        {{
            put("..", "\\.\\.");
            put("::", "::");
        }};
        String upper = symbol.toUpperCase();

        for (Map.Entry<String, String> separatorEntry : groupSeparators.entrySet())
        {
            Set<String> geneIDs = new HashSet<>();
            if (upper.contains(separatorEntry.getKey()))
            {
                String[] parts = upper.split(separatorEntry.getValue());
                for (String part : parts)
                {
                    geneIDs.addAll(symbolToEnsg(part));
                }
                return geneIDs;
            }
        }

        if (symbolEnsg.containsKey(upper))
        {
            return symbolEnsg.get(upper);
        } else if (upper.matches("[A-Z\\d]+\\.VAR\\.\\d+.*"))
        {
            logger.debug("Removed VAR part from " + symbol);
            return symbolEnsg.get(upper.substring(0, upper.indexOf(".")));
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

    public String ensgToDescription(String ensg) throws NoSuchFieldException
    {
        if (ensgDescription.containsKey(ensg.toUpperCase()))
        {
            return ensgDescription.get(ensg.toUpperCase());
        } else
        {
            throw new NoSuchFieldException("Could not find ensg " + ensg + " in the map file.");
        }
    }

    public String symbolToDescription(String symbol) throws NoSuchFieldException
    {
        if (symbolDescription.containsKey(symbol.toUpperCase()))
        {
            return symbolDescription.get(symbol.toUpperCase());
        } else
        {
            throw new NoSuchFieldException("Could not find symbol " + symbol + " in the map file.");
        }
    }

    public String anyToDescription(String any) throws NoSuchFieldException
    {
        if (symbolDescription.containsKey(any.toUpperCase()))
        {
            return symbolDescription.get(any.toUpperCase());
        } else if (ensgDescription.containsKey(any.toUpperCase()))
        {
            return ensgDescription.get(any.toUpperCase());
        } else
        {
            throw new NoSuchFieldException("Could not find symbol or ensg " + any + " in the map file.");
        }
    }

    public boolean hasEnsg(String ensg)
    {
        return ensgSymbol.containsKey(ensg.toUpperCase());
    }

    public boolean hasSymbol(String symbol)
    {
        return symbolEnsg.containsKey(symbol.toUpperCase());
    }

    public String getEnsgByIndex(int index)
    {
        return ensgList.get(index);
    }
}
