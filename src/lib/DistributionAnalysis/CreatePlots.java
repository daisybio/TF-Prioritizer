package lib.DistributionAnalysis;

import lib.ExecutableStep;
import tfprio.TFPRIO;
import util.Configs.Config;
import util.FileFilters.Filters;

import java.io.*;
import java.util.*;

import static util.FileManagement.*;
import static util.ScriptExecution.executeAndWait;

public class CreatePlots extends ExecutableStep
{
    private final Config<File> d_input_tfDistributionHm =
            TFPRIO.configs.distributionAnalysis.fileStructure.d_tfTgScores_tfDistribution_hm;
    private final Config<File> d_input_tfDistributionAll =
            TFPRIO.configs.distributionAnalysis.fileStructure.d_tfTgScores_tfDistribution_all;
    private final Config<File> d_input_backgroundDistributionHm =
            TFPRIO.configs.distributionAnalysis.fileStructure.d_tfTgScores_backgroundDistribution_hm;
    private final Config<File> d_input_backgroundDistributionAll =
            TFPRIO.configs.distributionAnalysis.fileStructure.d_tfTgScores_backgroundDistribution_all;
    private final Config<File> f_plotScriptTemplate = TFPRIO.configs.scriptTemplates.f_distributionPlots;
    private final Config<File> f_mwuScriptTemplate = TFPRIO.configs.scriptTemplates.f_distributionMwuPlots;

    private final Config<File> d_output_plotsHm = TFPRIO.configs.distributionAnalysis.fileStructure.d_plots_hm;
    private final Config<File> d_output_plotsAll = TFPRIO.configs.distributionAnalysis.fileStructure.d_plots_all;
    private final Config<File> d_output_statsHm = TFPRIO.configs.distributionAnalysis.fileStructure.d_stats_hm;
    private final Config<File> d_output_statsAll = TFPRIO.configs.distributionAnalysis.fileStructure.d_stats_all;
    private final Config<File> d_output_scriptsHm = TFPRIO.configs.distributionAnalysis.fileStructure.d_plotsScripts_hm;
    private final Config<File> d_output_scriptsAll =
            TFPRIO.configs.distributionAnalysis.fileStructure.d_plotsScripts_all;
    private final Config<File> d_output_mwuScriptsHm = TFPRIO.configs.distributionAnalysis.fileStructure.d_mwuScriptsHm;
    private final Config<File> d_output_mwuScriptsAll =
            TFPRIO.configs.distributionAnalysis.fileStructure.d_mwuScriptsAll;
    private final Config<File> d_output_mwuPlots = TFPRIO.configs.distributionAnalysis.fileStructure.d_mwuPlots;

    private final Config<String> s_plotScript =
            TFPRIO.configs.distributionAnalysis.fileStructure.s_plotsScripts_pythonScript;
    private final Config<String> s_mwuScript = TFPRIO.configs.distributionAnalysis.fileStructure.s_mwuScripts_script;
    private final Config<String> s_input =
            TFPRIO.configs.distributionAnalysis.fileStructure.s_tfTgScores_backgroundDistribution_csv;
    private final Config<String> s_stats = TFPRIO.configs.distributionAnalysis.fileStructure.s_stats_csv;

    @Override protected Set<Config<File>> getRequiredFileStructure()
    {
        return new HashSet<>();
    }

    @Override protected Set<Config<File>> getCreatedFileStructure()
    {
        return new HashSet<>();
    }

    @Override protected Set<Config<?>> getRequiredConfigs()
    {
        return new HashSet<>();
    }

    @Override protected void execute()
    {
        for (String hm : TFPRIO.existingHms)
        {
            File plotScript = createPlotScript(hm);
            //File mwuScript = createMwuScript(hm);
            executorService.submit(() -> executeAndWait(plotScript, logger));
            //executorService.submit(() -> executeAndWait(mwuScript, logger));

            // TODO: Implement generation of ALL plots
        }
    }

    private File createMwuScript(String hm)
    {
        String script = null;
        try
        {
            script = readFile(f_mwuScriptTemplate.get());
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }
        assert script != null;

        File f_background = extend(d_input_backgroundDistributionHm.get(), hm, s_input.get());
        File d_data = extend(d_input_tfDistributionHm.get(), hm);

        script = script.replace("{BACKGROUNDFILE}", f_background.getAbsolutePath());

        StringBuilder sb_calls = new StringBuilder();

        for (File f_data : Objects.requireNonNull(d_data.listFiles(Filters.fileFilter)))
        {
            String tfName = f_data.getName().substring(0, f_data.getName().lastIndexOf("."));
            sb_calls.append("generate('").append(f_data.getAbsolutePath()).append("', '").append(tfName).append("')\n");
        }

        script = script.replace("{ CALLS }", sb_calls.toString());

        File f_script = extend(d_output_mwuScriptsHm.get(), hm + ".R");

        try
        {
            writeFile(f_script, script);
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        return f_script;
    }

    private File createPlotScript(String hm)
    {
        String script = null;
        try
        {
            script = readFile(f_plotScriptTemplate.get());
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }
        assert script != null;

        File f_background = extend(d_input_backgroundDistributionHm.get(), hm, s_input.get());
        File d_data = extend(d_input_tfDistributionHm.get(), hm);
        File d_plots = extend(d_output_plotsHm.get(), hm);
        File f_stats = extend(d_output_statsHm.get(), hm, s_stats.get());

        makeSureFileExists(f_stats, logger);
        makeSureDirectoryExists(d_plots, logger);

        script = script.replace("{BACKGROUNDFILE}", f_background.getAbsolutePath());
        script = script.replace("{STATSFILE}", f_stats.getAbsolutePath());

        StringBuilder sb_calls = new StringBuilder();

        for (File f_data : Objects.requireNonNull(d_data.listFiles(Filters.fileFilter)))
        {
            String tfName = f_data.getName().substring(0, f_data.getName().lastIndexOf("."));

            File f_plot = extend(d_plots, tfName + ".png");

            sb_calls.append("generate('").append(f_data.getAbsolutePath()).append("', '")
                    .append(f_plot.getAbsolutePath()).append("', '").append(tfName).append("')\n");
        }

        script = script.replace("{CALLS}", sb_calls.toString());

        File f_script = extend(d_output_scriptsHm.get(), hm + ".py");

        try
        {
            writeFile(f_script, script);
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        return f_script;
    }
}