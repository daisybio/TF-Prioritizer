package lib.Plots;

import lib.ExecutableStep;
import tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import static util.FileManagement.*;
import static util.ScriptExecution.executeAndWait;

public class GroupPlots extends ExecutableStep
{
    private final AbstractConfig<File> d_input = TFPRIO.configs.dynamite.fileStructure.d_output;
    private final AbstractConfig<File> f_scriptTemplate = TFPRIO.configs.scriptTemplates.f_plots_groupPlots;

    private final AbstractConfig<File> d_output = TFPRIO.configs.plots.fileStructure.d_output;
    private final AbstractConfig<File> d_outputData = TFPRIO.configs.plots.fileStructure.d_data;

    private final AbstractConfig<String> s_input = TFPRIO.configs.dynamite.fileStructure.s_output_toBePlotted;
    private final AbstractConfig<String> s_allDataSame = TFPRIO.configs.plots.fileStructure.s_data_hmLevelSame;
    private final AbstractConfig<String> s_allDataDifferent =
            TFPRIO.configs.plots.fileStructure.s_data_hmLevelDifferent;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(Arrays.asList(d_input, f_scriptTemplate));
    }

    @Override public Set<AbstractConfig<File>> getCreatedFileStructure()
    {
        return new HashSet<>(Arrays.asList(d_output, d_outputData));
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>(Arrays.asList(s_input, s_allDataSame, s_allDataDifferent));
    }

    @Override protected void execute()
    {
        String scriptTemplate = null;
        try
        {
            scriptTemplate = readFile(f_scriptTemplate.get());
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        for (String hm : TFPRIO.existingHms)
        {
            for (double threshold : TFPRIO.configs.plots.thresholds.get())
            {
                StringBuilder sb_calls = new StringBuilder();

                HashSet<String> pairings_differentStages = new HashSet<>();
                HashSet<String> pairings_sameStages = new HashSet<>();

                for (String pairing : TFPRIO.groupCombinationsToHms.keySet())
                {
                    if (!TFPRIO.groupCombinationsToHms.get(pairing).contains(hm))
                    {
                        continue;
                    }

                    String[] group_split = pairing.split("_");
                    String group1 = group_split[0];
                    String group2 = group_split[1];
                    if (group1.charAt(0) != group2.charAt(0))
                    {
                        pairings_differentStages.add(pairing);
                    } else
                    {
                        pairings_sameStages.add(pairing);
                    }

                    File f_input = extend(d_input.get(), pairing, hm, s_input.get());
                    File f_output = extend(d_output.get(), hm, String.valueOf(threshold), pairing + ".png");

                    sb_calls.append(pairing).append("=generate('").append(f_input.getAbsolutePath()).append("', ")
                            .append(threshold).append(", '").append(f_output.getAbsolutePath()).append("', '")
                            .append(hm).append("', '").append(group1).append("', '").append(group2).append("')")
                            .append("\n");
                }

                String script = scriptTemplate.replace("{CALLS}", sb_calls.toString());
                File f_script = extend(d_output.get(), hm, String.valueOf(threshold), "script.py");
                File f_outputData_same = extend(d_outputData.get(), hm, String.valueOf(threshold), s_allDataSame.get());
                File f_outputData_different =
                        extend(d_outputData.get(), hm, String.valueOf(threshold), s_allDataDifferent.get());

                File f_plotSame = extend(d_output.get(), hm, String.valueOf(threshold), "same_stages.png");
                File f_plotDifferent = extend(d_output.get(), hm, String.valueOf(threshold), "different_stages.png");

                makeSureFileExists(f_outputData_same, logger);

                script = script.replace("{DIFFERENT_STAGES}",
                        getStageCommands("different", f_outputData_different, f_plotDifferent,
                                pairings_differentStages));

                script = script.replace("{SAME_STAGES}",
                        getStageCommands("same", f_outputData_same, f_plotSame, pairings_sameStages));

                try
                {
                    writeFile(f_script, script);
                } catch (IOException e)
                {
                    logger.error(e.getMessage());
                }

                String command = "python3 " + f_script.getAbsolutePath();


                executorService.submit(() -> executeAndWait(command, logger));
            }
        }
    }

    private String getStageCommands(String name, File f_data, File f_plot, Set<String> pairings)
    {
        StringBuilder sb_commands = new StringBuilder();

        for (String pairing : pairings)
        {
            sb_commands.append("if('Peak' in ");
            sb_commands.append(pairing);
            sb_commands.append(".index):\n");
            sb_commands.append(" ".repeat(4));
            sb_commands.append(pairing);
            sb_commands.append("=pd.DataFrame.drop(");
            sb_commands.append(pairing);
            sb_commands.append(",index='Peak')\n");

            sb_commands.append(pairing);
            sb_commands.append("=");
            sb_commands.append(pairing);
            sb_commands.append("[~").append(pairing).append(".index.duplicated(keep='first')]\n");
        }

        sb_commands.append(name).append("=pd.concat([").append(String.join(", ", pairings)).append("], axis=1)\n");

        sb_commands.append("stages('").append(f_data.getAbsolutePath()).append("', '").append(f_plot.getAbsolutePath())
                .append("', ").append(name).append(")\n");

        return sb_commands.toString();
    }
}
