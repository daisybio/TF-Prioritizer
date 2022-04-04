package lib.DistributionAnalysis;

import lib.ExecutableStep;
import tfprio.TFPRIO;
import util.Configs.Config;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static util.FileManagement.readFile;
import static util.FileManagement.writeFile;
import static util.ScriptExecution.executeAndWait;

public class GenerateTargetGenesHeatmaps extends ExecutableStep
{
    private final Config<File> f_scriptTemplate = TFPRIO.configs.scriptTemplates.f_heatmaps;
    private final Config<File> f_map = TFPRIO.configs.deSeq2.fileStructure.f_mapping;
    private final Config<File> d_deseq2Preprocessing = TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_combined;
    private final Config<File> d_targetGenes = TFPRIO.configs.distributionAnalysis.fileStructure.d_dcg_targetGenes;

    private final Config<File> d_output = TFPRIO.configs.distributionAnalysis.fileStructure.d_heatmaps;
    private final Config<File> f_script = TFPRIO.configs.distributionAnalysis.fileStructure.f_heatmaps_script;

    private final Config<Integer> kTargetGenes = TFPRIO.configs.plots.topKGenes;

    @Override protected Set<Config<File>> getRequiredFileStructure()
    {
        return new HashSet<>(Arrays.asList(f_scriptTemplate, f_map, d_deseq2Preprocessing, d_targetGenes));
    }

    @Override protected Set<Config<File>> getCreatedFileStructure()
    {
        return new HashSet<>(Arrays.asList(d_output, f_script));
    }

    @Override protected Set<Config<?>> getRequiredConfigs()
    {
        return new HashSet<>(List.of(kTargetGenes));
    }

    @Override protected void execute()
    {
        String script = null;
        try
        {
            script = readFile(f_scriptTemplate.get());
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }
        assert script != null;

        script = script.replace("{TARGET_GENES_PATH}", d_targetGenes.get().getAbsolutePath());
        script = script.replace("{MAP_PATH}", f_map.get().getAbsolutePath());
        script = script.replace("{PREPROCESSING_PATH}", d_deseq2Preprocessing.get().getAbsolutePath());
        script = script.replace("{HEATMAP_DIR}", d_output.get().getAbsolutePath());

        script = script.replace("{ NUMBER_OF_GENES }", String.valueOf(kTargetGenes.get()));

        try
        {
            writeFile(f_script.get(), script);
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        executeAndWait(f_script.get(), logger);
    }
}
