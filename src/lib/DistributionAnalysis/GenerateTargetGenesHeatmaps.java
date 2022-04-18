package lib.DistributionAnalysis;

import lib.ExecutableStep;
import tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;

import java.io.File;
import java.io.IOException;
import java.util.*;

import static util.FileManagement.readFile;
import static util.FileManagement.writeFile;
import static util.ScriptExecution.executeAndWait;

public class GenerateTargetGenesHeatmaps extends ExecutableStep
{
    private final AbstractConfig<File> f_scriptTemplate = TFPRIO.configs.scriptTemplates.f_heatmaps;
    private final AbstractConfig<File> f_map = TFPRIO.configs.deSeq2.fileStructure.f_mapping;
    private final AbstractConfig<File> d_deseq2Preprocessing =
            TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_combined;
    private final AbstractConfig<File> d_targetGenes =
            TFPRIO.configs.distributionAnalysis.fileStructure.d_dcg_targetGenes;

    private final GeneratedFileStructure d_output = TFPRIO.configs.distributionAnalysis.fileStructure.d_heatmaps;
    private final GeneratedFileStructure f_script = TFPRIO.configs.distributionAnalysis.fileStructure.f_heatmaps_script;

    private final AbstractConfig<Integer> kTargetGenes = TFPRIO.configs.plots.topKGenes;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(Arrays.asList(f_scriptTemplate, f_map, d_deseq2Preprocessing, d_targetGenes));
    }

    @Override public Set<GeneratedFileStructure> getCreatedFileStructure()
    {
        return new HashSet<>(Arrays.asList(d_output, f_script));
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
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

        StringBuilder sb_groups = new StringBuilder();
        StringBuilder sb_samples = new StringBuilder();
        StringBuilder sb_batches = new StringBuilder(", batch = c(");
        for (Map.Entry<String, String> entry : TFPRIO.sample_group.entrySet())
        {
            sb_samples.append("'").append(entry.getKey()).append("', ");
            sb_groups.append("'").append(entry.getValue()).append("', ");
            if (TFPRIO.sample_batch != null)
            {
                sb_batches.append("'").append(TFPRIO.sample_batch.get(entry.getKey())).append("', ");
            }
        }
        sb_groups.setLength(sb_groups.length() - 2);
        sb_samples.setLength(sb_samples.length() - 2);
        if (TFPRIO.sample_batch != null)
        {
            sb_batches.setLength(sb_batches.length() - 2);
            sb_batches.append(")");
            script = script.replace("{ BATCHES }", sb_batches.toString());
        } else
        {
            script = script.replace("{ BATCHES }", "");
        }
        script = script.replace("{ GROUPS }", sb_groups.toString());
        script = script.replace("{ SAMPLES }", sb_samples.toString());
        script = script.replace("{ DESIGN }", TFPRIO.sample_batch != null ? "batch+group" : "group");

        script = script.replace("{ NUMBER_OF_GENES }", String.valueOf(kTargetGenes.get()));

        try
        {
            writeFile(f_script.get(), script);
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        // TODO: Implement multithreading

        executeAndWait(f_script.get(), logger);
    }
}
