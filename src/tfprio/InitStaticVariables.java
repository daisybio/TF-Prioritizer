package tfprio;

import lib.ExecutableStep;
import util.Configs.Config;
import util.MapSymbolAndEnsg;

import java.io.File;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

public class InitStaticVariables extends ExecutableStep
{
    // Configs for map
    private final Config<File> f_map = TFPRIO.configs.deSeq2.fileStructure.f_mapping;
    private final Config<File> f_scriptTemplate = TFPRIO.configs.scriptTemplates.f_mapping;
    private final Config<File> f_script = TFPRIO.configs.deSeq2.fileStructure.f_mappingScript;
    private final Config<File> f_geneIDs = TFPRIO.configs.deSeq2.inputGeneID;
    private final Config<String> datasetSpecies = TFPRIO.configs.deSeq2.biomartDatasetSpecies;
    private final Config<String> datasetSymbolColumn = TFPRIO.configs.deSeq2.biomartDatasetSymbolColumn;

    // Configs for inputDirectory
    private final Config<File> inputDirectory = TFPRIO.configs.tepic.inputDirectory;

    @Override protected Set<Config<File>> getRequiredFileStructure()
    {
        return new HashSet<>(Arrays.asList(f_scriptTemplate, f_geneIDs, inputDirectory));
    }

    @Override protected Set<Config<File>> getCreatedFileStructure()
    {
        return new HashSet<>(Arrays.asList(f_map, f_script));
    }

    @Override protected Set<Config<?>> getRequiredConfigs()
    {
        return new HashSet<>(Arrays.asList(datasetSpecies, datasetSymbolColumn));
    }

    @Override protected void execute()
    {
        TFPRIO.latestInputDirectory = TFPRIO.configs.tepic.inputDirectory;
        TFPRIO.mapSymbolAndEnsg = new MapSymbolAndEnsg();
    }
}
