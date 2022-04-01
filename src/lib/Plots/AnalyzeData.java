package lib.Plots;

import lib.ExecutableStep;
import tfprio.TFPRIO;
import util.Configs.Config;
import util.FileFilters.Filters;

import java.io.*;
import java.util.*;

import static util.FileManagement.extend;
import static util.FileManagement.makeSureFileExists;

public class AnalyzeData extends ExecutableStep
{
    private final Config<File> d_input_readCounts = TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_meanCounts;
    private final Config<File> d_input_tpm = TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_tpm_tpmResults;
    private final Config<File> d_input_plotData = TFPRIO.configs.plots.fileStructure.d_data;
    private final Config<File> f_input_tfs = TFPRIO.configs.tepic.fileStructure.f_postprocessing_tfs_csv;

    private final Config<File> d_output_groupLevel = TFPRIO.configs.plots.fileStructure.d_analysisData_tpLevel;
    private final Config<File> d_output_hmLevel = TFPRIO.configs.plots.fileStructure.d_analysisData_hmLevel;
    private final Config<File> d_output_website = TFPRIO.configs.plots.fileStructure.d_analysisData_websiteOverview;

    private final Config<String> s_plotData_different = TFPRIO.configs.plots.fileStructure.s_data_hmLevelDifferent;
    private final Config<String> s_plotData_same = TFPRIO.configs.plots.fileStructure.s_data_hmLevelSame;

    @Override protected Set<Config<File>> getRequiredFileStructure()
    {
        return new HashSet<>(Arrays.asList(d_input_readCounts, d_input_tpm, d_input_plotData, f_input_tfs));
    }

    @Override protected Set<Config<File>> getCreatedFileStructure()
    {
        return new HashSet<>(Arrays.asList(d_output_groupLevel, d_output_hmLevel, d_output_website));
    }

    @Override protected Set<Config<?>> getRequiredConfigs()
    {
        return new HashSet<>(Arrays.asList(s_plotData_different, s_plotData_same));
    }

    @Override protected void execute()
    {
        HashMap<String, HashMap<String, Double>> group_geneID_tpm = new HashMap<>();

        for (File f_group : Objects.requireNonNull(d_input_tpm.get().listFiles(Filters.fileFilter)))
        {
            String group = f_group.getName().substring(0, f_group.getName().lastIndexOf("."));

            HashMap<String, Double> geneID_tpm = new HashMap<>();

            try (BufferedReader br = new BufferedReader(new FileReader(f_group)))
            {
                String line;
                br.readLine();

                while ((line = br.readLine()) != null)
                {
                    String[] split = line.split("\t");
                    geneID_tpm.put(split[1], Double.parseDouble(split[4]));
                }
            } catch (IOException e)
            {
                logger.error(e.getMessage());
            }

            group_geneID_tpm.put(group, geneID_tpm);
        }

        HashMap<Double, HashMap<String, HashMap<String, HashMap<String, Double>>>> cutoff_hm_tf_counts_DIFFERENT =
                new HashMap<>();
        HashMap<Double, HashMap<String, HashMap<String, HashMap<String, Double>>>> cutoff_hm_tf_counts_SAME =
                new HashMap<>();
        HashMap<String, String> composed_tfs = new HashMap<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(f_input_tfs.get())))
        {
            String inputLine;

            while ((inputLine = reader.readLine()) != null)
            {
                String[] split = inputLine.split("\t");

                for (int i = 1; i < split.length; i++)
                {
                    composed_tfs.put(split[i].toUpperCase(), split[0].toUpperCase());
                }
            }
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        HashMap<String, HashMap<String, Double>> group_tf_counts = new HashMap<>();

        for (File d_group : Objects.requireNonNull(d_input_readCounts.get().listFiles(Filters.fileFilter)))
        {
            String group = d_group.getName().substring(0, d_group.getName().lastIndexOf("."));

            HashMap<String, Double> tf_counts = new HashMap<>();

            try (BufferedReader br = new BufferedReader(new FileReader(d_group)))
            {
                String line;
                br.readLine();

                while ((line = br.readLine()) != null)
                {
                    String[] split = line.split("\t");
                    String geneSymbol = split[0].toUpperCase();

                    if (composed_tfs.containsKey(geneSymbol))
                    {
                        String composedSymbol = composed_tfs.get(geneSymbol);

                        if (!tf_counts.containsKey(composedSymbol))
                        {
                            tf_counts.put(composedSymbol, .0);
                        }

                        double count = Double.parseDouble(split[2]);
                        tf_counts.put(composedSymbol, tf_counts.get(composedSymbol) + count);
                    } else
                    {
                        tf_counts.put(split[0].toUpperCase(), Double.parseDouble(split[2]));
                    }
                }
            } catch (IOException e)
            {
                logger.error(e.getMessage());
            }

            group_tf_counts.put(group, tf_counts);
        }


        HashMap<String, HashMap<Double, HashMap<String, Boolean>>> hm_th_tf_found = new HashMap<>();

        for (String hm : TFPRIO.existingHms)
        {
            for (double threshold : TFPRIO.configs.plots.thresholds.get())
            {
                HashMap<Double, HashMap<String, Boolean>> th_tf_found;

                HashMap<String, Boolean> tf_found;

                if (hm_th_tf_found.containsKey(hm))
                {
                    th_tf_found = hm_th_tf_found.get(hm);
                } else
                {
                    th_tf_found = new HashMap<>();
                }

                if (th_tf_found.containsKey(threshold))
                {
                    tf_found = th_tf_found.get(threshold);
                } else
                {
                    tf_found = new HashMap<>();

                    for (Object prioTf : TFPRIO.configs.igv.importantLociAllPrioTf.get())
                    {
                        tf_found.put((String) prioTf, false);
                    }
                    th_tf_found.put(threshold, tf_found);
                }

                hm_th_tf_found.put(hm, th_tf_found);

                File f_availableTfs = extend(d_output_website.get(), hm, String.valueOf(threshold),
                        TFPRIO.configs.plots.fileStructure.s_analysisData_websiteOverview_tfAvailable.get());
                File f_different =
                        extend(d_input_plotData.get(), hm, String.valueOf(threshold), s_plotData_different.get());
                File f_same = extend(d_input_plotData.get(), hm, String.valueOf(threshold), s_plotData_same.get());

                makeSureFileExists(f_availableTfs, logger);
                makeSureFileExists(f_different, logger);
                makeSureFileExists(f_same, logger);

                processPlotData(f_different, threshold, hm, tf_found, group_tf_counts, group_geneID_tpm,
                        cutoff_hm_tf_counts_DIFFERENT);
                processPlotData(f_same, threshold, hm, tf_found, group_tf_counts, group_geneID_tpm,
                        cutoff_hm_tf_counts_SAME);

                try (BufferedWriter writer = new BufferedWriter(new FileWriter(f_availableTfs)))
                {
                    writer.write("TF\tAVAILABLE");
                    writer.newLine();
                    for (String tf_key : tf_found.keySet())
                    {
                        writer.write(tf_key + "\t" + tf_found.get(tf_key));
                        writer.newLine();
                    }
                } catch (IOException e)
                {
                    logger.error(e.getMessage());
                }
            }
        }

        for (Double threshold : cutoff_hm_tf_counts_DIFFERENT.keySet())
        {
            File f_hmLevel = extend(d_output_hmLevel.get(), String.valueOf(threshold), s_plotData_different.get());
            makeSureFileExists(f_hmLevel, logger);

            HashMap<String, HashMap<String, HashMap<String, Double>>> available_hms =
                    cutoff_hm_tf_counts_DIFFERENT.get(threshold);

            ArrayList<HashMap<String, HashMap<String, Double>>> tf_lists = new ArrayList<>();

            for (String kk : available_hms.keySet())
            {
                tf_lists.add(available_hms.get(kk));
            }

            HashMap<String, HashMap<String, Double>> results = new HashMap<>();

            for (int i = 0; i < tf_lists.size(); i++)
            {
                for (String tf : tf_lists.get(i).keySet())
                {
                    int count_hms = 0;

                    for (HashMap<String, HashMap<String, Double>> list : tf_lists)
                    {
                        if (list.containsKey(tf))
                        {
                            count_hms++;
                        }
                    }

                    if (count_hms >= TFPRIO.configs.plots.cutoffHms.get())
                    {
                        results.put(tf, tf_lists.get(i).get(tf));
                    }
                }
            }


            try (BufferedWriter bw = new BufferedWriter(new FileWriter(f_hmLevel)))
            {
                StringBuilder sb = new StringBuilder();
                ArrayList<String> tp_order = new ArrayList<>();
                sb.append("TF");
                for (String tp : group_tf_counts.keySet())
                {
                    sb.append("\t");
                    sb.append(tp);
                    tp_order.add(tp);
                }
                bw.write(sb.toString());
                bw.newLine();

                for (String tf : results.keySet())
                {
                    StringBuilder sb_tf = new StringBuilder();
                    sb_tf.append(tf);
                    HashMap<String, Double> tf_info = results.get(tf);

                    for (String order : tp_order)
                    {
                        sb_tf.append("\t");
                        sb_tf.append(tf_info.get(order));
                    }
                    bw.write(sb_tf.toString());
                    bw.newLine();
                }
            } catch (IOException e)
            {
                logger.error(e.getMessage());
            }
        }
    }

    private void processPlotData(File source, double threshold, String hm, HashMap<String, Boolean> tf_found,
                                 HashMap<String, HashMap<String, Double>> group_tf_counts,
                                 HashMap<String, HashMap<String, Double>> group_geneID_tpm,
                                 HashMap<Double, HashMap<String, HashMap<String, HashMap<String, Double>>>> data)
    {
        logger.debug("Processing plot data: " + source.getAbsolutePath());

        if (!data.containsKey(threshold))
        {
            data.put(threshold, new HashMap<>());
        }

        if (!data.get(threshold).containsKey(hm))
        {
            data.get(threshold).put(hm, new HashMap<>());
        }

        File targetFile = extend(d_output_groupLevel.get(), hm, String.valueOf(threshold), source.getName());
        makeSureFileExists(targetFile, logger);
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(targetFile));
             BufferedReader reader = new BufferedReader(new FileReader(source)))
        {
            StringBuilder sb = new StringBuilder();
            sb.append("TF");
            for (String group : group_tf_counts.keySet())
            {
                sb.append("\t");
                sb.append(group);
            }
            writer.write(sb.toString());
            writer.newLine();


            String line;
            reader.readLine();
            while ((line = reader.readLine()) != null)
            {
                String[] split = line.split(",");

                if (tf_found.containsKey(split[0].toUpperCase()))
                {
                    tf_found.put(split[0].toUpperCase(), true);
                }

                int count = 0;
                for (String s : split)
                {
                    if (!s.equals(""))
                    {
                        count++;
                    }

                }

                // TODO: Handling of non-mappable genes has to be improved
                Set<String> geneIDs;
                try
                {
                    geneIDs = TFPRIO.mapSymbolAndEnsg.symbolToEnsg(split[0]);
                } catch (NoSuchFieldException e)
                {
                    logger.warn(e.getMessage());
                    geneIDs = null;
                }

                assert geneIDs != null;

                if (count > TFPRIO.configs.plots.cutoffTps.get())
                {
                    HashMap<String, Double> tf_gc = new HashMap<>();
                    //tf_gc.put(split[0],group_tf_counts.get(name).get(split[0]));
                    StringBuilder sb_intern = new StringBuilder();
                    sb_intern.append(split[0]);
                    int count_row = 0;

                    for (String k : group_tf_counts.keySet())
                    {
                        if (group_tf_counts.get(k).containsKey(split[0].toUpperCase()))
                        {
                            count_row += group_tf_counts.get(k).get(split[0].toUpperCase());

                            tf_gc.put(k, group_tf_counts.get(k).get(split[0].toUpperCase()));
                            sb_intern.append("\t");
                            sb_intern.append(group_tf_counts.get(k).get(split[0].toUpperCase()));
                        }
                    }

                    double cumm_tpm = 0.0;

                    for (String key_group : group_geneID_tpm.keySet())
                    {
                        HashMap<String, Double> lookup = group_geneID_tpm.get(key_group);
                        for (String geneID : geneIDs)
                        {
                            if (lookup.containsKey(geneID))
                            {
                                cumm_tpm += lookup.get(geneID);
                                break;
                            }
                        }
                    }

                    boolean passed_count_threshold = count_row >= TFPRIO.configs.plots.cutoffGcs.get();
                    boolean passed_tpm_threshold = cumm_tpm >= TFPRIO.configs.plots.cutoffTpms.get();

                    if (passed_count_threshold && passed_tpm_threshold)
                    {
                        data.get(threshold).get(hm).put(split[0].toUpperCase(), tf_gc);
                        writer.write(sb_intern.toString());
                        writer.newLine();
                    }
                }
            }
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }
    }
}
