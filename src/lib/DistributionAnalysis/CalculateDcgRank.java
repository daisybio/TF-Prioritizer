package lib.DistributionAnalysis;

import lib.DistributionAnalysis.Classes.Stats;
import lib.DistributionAnalysis.Classes.StatsCollection;
import lib.DistributionAnalysis.Classes.StatsCummulativeGain;
import lib.ExecutableStep;
import tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.FileFilters.Filters;

import java.io.*;
import java.util.*;

import static util.FileManagement.extend;
import static util.FileManagement.makeSureFileExists;

public class CalculateDcgRank extends ExecutableStep
{
    private final AbstractConfig<File> d_input_statsAll = TFPRIO.configs.distributionAnalysis.fileStructure.d_stats_all;
    private final AbstractConfig<File> d_input_statsHm = TFPRIO.configs.distributionAnalysis.fileStructure.d_stats_hm;

    private final AbstractConfig<File> f_output_dcg = TFPRIO.configs.distributionAnalysis.fileStructure.f_dcg_stats;

    private final AbstractConfig<Boolean> performAllAnalysis = TFPRIO.configs.distributionAnalysis.performAllAnalysis;
    private final AbstractConfig<String> s_stats = TFPRIO.configs.distributionAnalysis.fileStructure.s_stats_csv;

    private final AbstractConfig<String> allName = TFPRIO.configs.distributionAnalysis.allName;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(List.of(d_input_statsHm))
        {{
            if (performAllAnalysis.get())
            {
                add(d_input_statsAll);
            }
        }};
    }

    @Override public Set<AbstractConfig<File>> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(f_output_dcg));
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>(Arrays.asList(performAllAnalysis, s_stats))
        {{
            if (performAllAnalysis.get())
            {
                add(allName);
            }
        }};
    }

    @Override protected void execute()
    {
        HashMap<String, HashMap<String, Integer>> hm_tf_rank = new HashMap<>();
        HashMap<String, StatsCollection> hmStats = new HashMap<>();

        for (File d_hm : Objects.requireNonNull(d_input_statsHm.get().listFiles(Filters.directoryFilter)))
        {
            File f_hm = extend(d_hm, s_stats.get());
            logger.debug(f_hm.getAbsolutePath());
            StatsCollection hmCollection = getStats(f_hm);
            hmStats.put(d_hm.getName(), hmCollection);
        }

        if (performAllAnalysis.get())
        {
            // TODO: Currently the all directory is used, not the according file. Will lead to an error.
            StatsCollection allStats = getStats(d_input_statsAll.get());

            hmStats.put(allName.get(), allStats);
        }

        for (String hm : hmStats.keySet())
        {
            HashMap<String, Integer> tf_to_ranks = new HashMap<>();
            int rank = 1;

            ArrayList<Stats> tf_rank = hmStats.get(hm).getTfs();

            for (Stats tf : tf_rank)
            {
                tf_to_ranks.put(tf.getLabel(), rank);
                rank++;
            }
            hm_tf_rank.put(hm, tf_to_ranks);
        }


        ArrayList<StatsCummulativeGain> cummulativeGains = new ArrayList<>();
        HashSet<String> already_calculated_tfs = new HashSet<>();

        for (String hm : hm_tf_rank.keySet())
        {
            HashMap<String, Integer> tf_rank = hm_tf_rank.get(hm);
            for (String tfSymbol : tf_rank.keySet())
            {
                if (already_calculated_tfs.contains(tfSymbol))
                {
                    continue;
                }

                double score = 0;

                double tfScore = (double) (tf_rank.size() - tf_rank.get(tfSymbol)) / tf_rank.size();
                score += tfScore;

                HashMap<String, Stats> tfStats = new HashMap<>();
                HashMap<String, Stats> tfBackground = new HashMap<>();

                StatsCollection currentHmStats = hmStats.get(hm);
                tfStats.put(hm, currentHmStats.getTfs().get(tf_rank.get(tfSymbol) - 1));
                tfBackground.put(hm, currentHmStats.getBackground());


                for (String hm2 : hm_tf_rank.keySet())
                {
                    if (hm2.equals(hm))
                    {
                        continue;
                    }
                    HashMap<String, Integer> tf_rank2 = hm_tf_rank.get(hm2);

                    if (tf_rank2.containsKey(tfSymbol))
                    {
                        StatsCollection currentHmStats2 = hmStats.get(hm2);
                        tfStats.put(hm2, currentHmStats2.getTfs().get(tf_rank2.get(tfSymbol) - 1));
                        tfBackground.put(hm2, currentHmStats2.getBackground());

                        double score2 = (double) (tf_rank2.size() - tf_rank2.get(tfSymbol)) / tf_rank2.size();
                        score += score2;
                    }
                }
                StatsCummulativeGain ac = new StatsCummulativeGain(tfSymbol, score, tfStats, tfBackground);

                cummulativeGains.add(ac);
                already_calculated_tfs.add(tfSymbol);
            }
        }

        Collections.sort(cummulativeGains);

        makeSureFileExists(f_output_dcg.get(), logger);

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(f_output_dcg.get())))
        {
            writer.write("RANK\tTF\tSCORE");
            writer.newLine();
            int rankIndex = 1;
            for (StatsCummulativeGain ac : cummulativeGains)
            {
                writer.write(String.valueOf(rankIndex));
                writer.write('\t');
                writer.write(ac.toString());
                writer.newLine();
                rankIndex++;
            }
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }
    }

    private StatsCollection getStats(File f_stats)
    {
        StatsCollection allConsideredTfs = null;

        try (BufferedReader reader = new BufferedReader(new FileReader(f_stats)))
        {
            String line;
            reader.readLine();
            line = reader.readLine();

            Stats background = new Stats(line);

            allConsideredTfs = new StatsCollection(background);

            while ((line = reader.readLine()) != null)
            {
                Stats tf = new Stats(line);
                allConsideredTfs.add(tf);
            }
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        assert allConsideredTfs != null;
        allConsideredTfs.sort();

        return allConsideredTfs;
    }
}
