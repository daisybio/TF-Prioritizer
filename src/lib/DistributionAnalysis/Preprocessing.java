package lib.DistributionAnalysis;

import lib.ExecutableStep;
import tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;

import java.io.*;
import java.util.*;

import static util.FileManagement.extend;
import static util.FileManagement.makeSureFileExists;

public class Preprocessing extends ExecutableStep
{
    private final AbstractConfig<File> d_input = TFPRIO.configs.plots.fileStructure.d_data;

    private final AbstractConfig<File> f_output = TFPRIO.configs.distributionAnalysis.fileStructure.f_analyzedTfs;

    private final AbstractConfig<List<Double>> thresholds = TFPRIO.configs.plots.thresholds;
    private final AbstractConfig<String> s_plotData_same = TFPRIO.configs.plots.fileStructure.s_data_hmLevelSame;
    private final AbstractConfig<String> s_plotData_different =
            TFPRIO.configs.plots.fileStructure.s_data_hmLevelDifferent;
    private final AbstractConfig<String> sameGroups = TFPRIO.configs.general.sameTps;
    private final AbstractConfig<String> differentGroups = TFPRIO.configs.general.differentTps;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(List.of(d_input));
    }

    @Override public Set<AbstractConfig<File>> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(f_output));
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>(
                Arrays.asList(thresholds, s_plotData_same, s_plotData_different, sameGroups, differentGroups));
    }

    @Override protected void execute()
    {
        HashMap<String, HashMap<String, HashSet<String>>> distinct_tf_hm_diff_same = new HashMap<>();

        Map<String, String> stage_name = new HashMap<>()
        {{
            put(s_plotData_different.get(), differentGroups.get());
            put(s_plotData_same.get(), sameGroups.get());
        }};

        for (String hm : TFPRIO.existingHms)
        {
            for (String suffix : stage_name.keySet())
            {
                File f_input = extend(d_input.get(), hm, String.valueOf(thresholds.get().get(0)), suffix);
                try (BufferedReader reader = new BufferedReader(new FileReader(f_input)))
                {
                    String line_stages;
                    reader.readLine();
                    while ((line_stages = reader.readLine()) != null)
                    {
                        String[] split = line_stages.split(",");
                        String tf = split[0].toUpperCase();

                        if (!distinct_tf_hm_diff_same.containsKey(tf))
                        {
                            distinct_tf_hm_diff_same.put(tf, new HashMap<>());
                        }

                        if (!distinct_tf_hm_diff_same.get(tf).containsKey(hm))
                        {
                            distinct_tf_hm_diff_same.get(tf).put(hm, new HashSet<>());
                        }
                        distinct_tf_hm_diff_same.get(tf).get(hm).add(stage_name.get(suffix));
                    }
                } catch (IOException e)
                {
                    logger.error(e.getMessage());
                }
            }
        }


        makeSureFileExists(f_output.get(), logger);
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(f_output.get())))
        {
            writer.write("TF\tHM\tSTAGES");
            writer.newLine();

            for (String tf : distinct_tf_hm_diff_same.keySet())
            {
                HashMap<String, HashSet<String>> hm_diff_stages = distinct_tf_hm_diff_same.get(tf);

                for (String hm : hm_diff_stages.keySet())
                {
                    writer.write(tf);
                    writer.write("\t");
                    writer.write(hm);
                    writer.write("\t");
                    writer.write(String.join(";", hm_diff_stages.get(hm)));
                    writer.newLine();
                }
            }
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }
    }
}
