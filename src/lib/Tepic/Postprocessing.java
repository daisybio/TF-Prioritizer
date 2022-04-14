package lib.Tepic;

import lib.ExecutableStep;
import tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.FileFilters.Filters;

import java.io.*;
import java.util.*;

import static util.FileFilters.Filters.getSuffixFilter;
import static util.FileManagement.*;

public class Postprocessing extends ExecutableStep
{
    private AbstractConfig<File> d_input;

    private final AbstractConfig<File> d_postprocessingInput =
            TFPRIO.configs.tepic.fileStructure.d_postprocessing_input;
    private final AbstractConfig<File> d_output = TFPRIO.configs.tepic.fileStructure.d_postprocessing_output;
    private final AbstractConfig<File> f_output_tfs = TFPRIO.configs.tepic.fileStructure.f_postprocessing_tfs_csv;

    private final AbstractConfig<Boolean> mutuallyExclusive = TFPRIO.configs.mixOptions.mutuallyExclusive;
    private final AbstractConfig<Boolean> originalDecay = TFPRIO.configs.tepic.originalDecay;
    private final AbstractConfig<Integer> tpmCutoff = TFPRIO.configs.tepic.tpmCutoff;
    private final AbstractConfig<Boolean> doNotGenerate = TFPRIO.configs.tepic.doNotGenerate;
    private final AbstractConfig<String> s_meanAffinities =
            TFPRIO.configs.tepic.fileStructure.s_postprocessing_output_meanAffinitiesDir;
    private final AbstractConfig<String> s_ratios =
            TFPRIO.configs.tepic.fileStructure.s_postprocessing_output_ratiosDir;


    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(List.of(d_input));
    }

    @Override protected Set<AbstractConfig<File>> getCreatedFileStructure()
    {
        return new HashSet<>(Arrays.asList(d_output, d_postprocessingInput, f_output_tfs));
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>(
                Arrays.asList(tpmCutoff, mutuallyExclusive, originalDecay, doNotGenerate, s_meanAffinities, s_ratios));
    }

    @Override protected void updateInputDirectory()
    {
        d_input = TFPRIO.configs.tepic.randomizeTfGeneMatrix.get() ?
                TFPRIO.configs.tepic.fileStructure.d_outputRaw_shuffle : TFPRIO.configs.tepic.fileStructure.d_outputRaw;
    }

    @Override protected int getShutDownTimeOutMinutes()
    {
        return 30;
    }

    @Override protected void execute()
    {
        Set<String> tf_groups = copyInputAndGetTfGroups();

        makeSureFileExists(f_output_tfs.get(), logger);

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(f_output_tfs.get())))
        {
            for (String tf_group : tf_groups)
            {
                String[] split = tf_group.split("[.][.]");

                writer.write(tf_group);

                for (String s : split)
                {
                    writer.write("\t");
                    writer.write(s.toUpperCase());
                }
                writer.newLine();
            }
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        //generate output structure
        //HMs

        for (Map.Entry<String, Set<String>> entry : TFPRIO.groupCombinationsToHms.entrySet())
        {
            for (String hm : entry.getValue())
            {
                File directory = extend(d_output.get(), entry.getKey(), hm);

                try
                {
                    makeSureDirectoryExists(directory);
                } catch (IOException e)
                {
                    logger.error(e.getMessage());
                }
            }
        }

        //prepare command line for all (computeMeanRatioTFAffinities.py
        HashSet<String> all_tfs = new HashSet<>();

        for (File d_pairing : Objects.requireNonNull(d_output.get().listFiles(Filters.directoryFilter)))
        {
            String pairing = d_pairing.getName();
            String[] pairingSplit = pairing.split("_");
            assert pairingSplit.length == 2;

            String group1 = pairingSplit[0];
            String group2 = pairingSplit[1];

            for (File d_hm : Objects.requireNonNull(d_pairing.listFiles(Filters.directoryFilter)))
            {
                String hm = d_hm.getName();

                File d_input1;
                File d_input2;

                if (!mutuallyExclusive.get())
                {
                    d_input1 = extend(d_postprocessingInput.get(), group1, hm);
                    d_input2 = extend(d_postprocessingInput.get(), group2, hm);
                } else
                {
                    d_input1 = extend(d_postprocessingInput.get(), pairing, hm);
                    d_input2 = extend(d_postprocessingInput.get(), pairing, hm);
                }

                File[] samples_group1_check_tfs = d_input1.listFiles();
                File[] samples_group2_check_tfs = d_input2.listFiles();
                assert samples_group1_check_tfs != null;
                assert samples_group2_check_tfs != null;


                try (BufferedReader br_group1_check_tfs = new BufferedReader(
                        new FileReader((samples_group1_check_tfs[0])));
                     BufferedReader br_group2_check_tfs = new BufferedReader(
                             new FileReader((samples_group2_check_tfs[0]))))
                {
                    String header_group1_check_tfs = br_group1_check_tfs.readLine();
                    all_tfs.addAll(Arrays.asList(header_group1_check_tfs.split("\t")));

                    String header_group2_check_tfs = br_group2_check_tfs.readLine();
                    all_tfs.addAll(Arrays.asList(header_group2_check_tfs.split("\t")));
                } catch (IOException e)
                {
                    logger.error(e.getMessage());
                }

                File group1_input_dir;
                File group2_input_dir;

                //if TPM filter was used we need to postprocess the TPM files. Otherwise, it will not work!
                if (tpmCutoff.get() > 0)
                {
                    logger.info("TPM filter > 0, start postprocessing of TPM filtered scores");
                    File d_output1 = extend(d_output.get(), hm, pairing, group1);
                    File d_output2 = extend(d_output.get(), hm, pairing, group2);

                    //build intersect and write new files with filter

                    File folder_group1;
                    File folder_group2;
                    File[] samples_group1;
                    File[] samples_group2;

                    if (!mutuallyExclusive.get())
                    {
                        folder_group1 = extend(d_postprocessingInput.get(), group1, hm);
                        folder_group2 = extend(d_postprocessingInput.get(), group2, hm);

                        samples_group1 = folder_group1.listFiles();
                        samples_group2 = folder_group2.listFiles();
                    } else
                    {
                        folder_group1 = extend(d_postprocessingInput.get(), pairing, hm);
                        folder_group2 = extend(d_postprocessingInput.get(), pairing, hm);
                        samples_group1 = new File[1];
                        samples_group1[0] = Objects.requireNonNull(folder_group1.listFiles())[0];
                        samples_group2 = new File[1];
                        samples_group2[0] = Objects.requireNonNull(folder_group2.listFiles())[1];
                    }

                    assert samples_group1 != null;
                    assert samples_group2 != null;

                    for (File f : Objects.requireNonNull(folder_group1.listFiles()))
                    {
                        softLink(extend(d_output1, f.getName()), f, logger);

                        if (mutuallyExclusive.get())
                        {
                            break;
                        }
                    }


                    for (File f : Objects.requireNonNull(folder_group2.listFiles()))
                    {
                        softLink(extend(d_output2, f.getName()), f, logger);

                        if (mutuallyExclusive.get())
                        {
                            break;
                        }
                    }
                    logger.info("TPM filter > 0, end postprocessing of TPM filtered scores");

                    // Set to post processed TMP filtered data
                    group1_input_dir = d_output1;
                    group2_input_dir = d_output2;
                } else
                {
                    if (!mutuallyExclusive.get())
                    {
                        group1_input_dir = extend(d_postprocessingInput.get(), group1, hm);
                        group2_input_dir = extend(d_postprocessingInput.get(), group2, hm);
                    } else
                    {
                        //NOTE: THIS IS NOT TESTED AND COULD CAUSE TROUBLE WITH NO TPM + MUTUALLY EXCLUSIVE OPTION
                        if (mutuallyExclusive.get())
                        {
                            File f_group1_data = Objects.requireNonNull(
                                    extend(d_postprocessingInput.get(), group1, hm).listFiles())[0];
                            File f_group2_data = Objects.requireNonNull(
                                    extend(d_postprocessingInput.get(), group1, hm).listFiles())[1];

                            File f_group1_file =
                                    extend(d_postprocessingInput.get(), pairing, hm, group1, f_group1_data.getName());
                            File f_group2_file =
                                    extend(d_postprocessingInput.get(), pairing, hm, group2, f_group2_data.getName());

                            softLink(f_group1_file, f_group1_data, logger);
                            softLink(f_group2_file, f_group2_data, logger);

                            group1_input_dir = f_group1_file;
                            group2_input_dir = f_group2_file;
                        } else
                        {
                            group1_input_dir = extend(d_postprocessingInput.get(), group1, hm);
                            group2_input_dir = extend(d_postprocessingInput.get(), group2, hm);
                        }
                    }
                }


                File output_mean_affinities = extend(d_output.get(), pairing, hm, s_meanAffinities.get());
                File output_ratios = extend(d_output.get(), pairing, hm, s_ratios.get());

                File f_output_meanAffinities1 = extend(output_mean_affinities, group1 + ".txt");
                File f_output_meanAffinities2 = extend(output_mean_affinities, group2 + ".txt");
                File f_output_ratios = extend(output_ratios, group1 + "_" + group2 + ".txt");

                executorService.submit(() -> calculate(group1_input_dir, group2_input_dir, f_output_meanAffinities1,
                        f_output_meanAffinities2, f_output_ratios));
            }
        }

        try (BufferedWriter bw_all_tfs = new BufferedWriter(new FileWriter(f_output_tfs.get(), true)))
        {
            for (String k : all_tfs)
            {
                bw_all_tfs.write(k);
                bw_all_tfs.newLine();
            }
        } catch (IOException e)
        {
            e.printStackTrace();
        }
    }

    private Set<String> copyInputAndGetTfGroups()
    {
        String suffix;

        if (tpmCutoff.get() > 0)
        {
            suffix = "_Gene_View_Filtered_TPM.txt";
        } else
        {
            suffix = "_Gene_View_Filtered.txt";
        }

        Set<String> tf_groups = new HashSet<>();
        for (File d_group : Objects.requireNonNull(d_input.get().listFiles(Filters.directoryFilter)))
        {
            for (File d_hm : Objects.requireNonNull(d_group.listFiles(Filters.directoryFilter)))
            {
                for (File d_sample : Objects.requireNonNull(d_hm.listFiles(Filters.directoryFilter)))
                {
                    for (File f_data : Objects.requireNonNull(d_sample.listFiles(Filters.fileFilter)))
                    {
                        String name = f_data.getName();
                        if (name.matches(".*" + suffix))
                        {
                            //check tfs
                            try (BufferedReader reader = new BufferedReader(new FileReader(f_data)))
                            {
                                String line = reader.readLine();
                                String[] split = line.split("\t");
                                for (String st : split)
                                {
                                    for (String sep : Arrays.asList("::", ".."))
                                    {
                                        if (st.contains(sep))
                                        {
                                            String tfGroupString = String.join("..", st.split(sep));
                                            tf_groups.add(tfGroupString);
                                        }
                                    }
                                }
                                File f_link = extend(d_postprocessingInput.get(), d_group.getName(), d_hm.getName(),
                                        f_data.getName());
                                softLink(f_link, f_data);
                            } catch (IOException e)
                            {
                                e.printStackTrace();
                            }
                        }
                    }
                }
            }
        }
        return tf_groups;
    }

    private String getFileSuffix()
    {
        if (doNotGenerate.get())
        {
            if (originalDecay.get())
            {
                if (tpmCutoff.get() > 0)
                {
                    return "Three_Peak_Based_Features_Affinity_Gene_View_Filtered_TPM.txt";
                } else
                {
                    return "Three_Peak_Based_Features_Affinity_Gene_View_Filtered.txt";
                }
            } else
            {
                if (tpmCutoff.get() > 0)
                {
                    return "Peak_Features_Affinity_Gene_View_Filtered_TPM.txt";
                } else
                {
                    return "Peak_Features_Affinity_Gene_View_Filtered.txt";
                }
            }
        } else
        {
            if (originalDecay.get())
            {
                if (tpmCutoff.get() > 0)
                {
                    return "Peak_Features_Affinity_Gene_View_Filtered_TPM.txt";
                } else
                {
                    return "Peak_Features_Affinity_Gene_View_Filtered.txt";
                }
            } else
            {
                if (tpmCutoff.get() > 0)
                {
                    return "Affinity_Gene_View_Filtered_TPM.txt";
                } else
                {
                    return "Affinity_Gene_View_Filtered.txt";
                }
            }
        }
    }

    private synchronized void mergeAffinityFile(File file, Map<String, Map<String, Double>> affinityMap,
                                                String headerTemplate)
    {
        logger.debug("Loading affinity file: " + file.getAbsolutePath());
        try (BufferedReader reader = new BufferedReader(new FileReader(file)))
        {
            String header = reader.readLine();

            if (!header.equals(headerTemplate))
            {
                logger.error("File header does not match template header.");
            }

            String[] headerSplit = headerTemplate.split("\t");

            String inputLine;
            while ((inputLine = reader.readLine()) != null)
            {
                String[] inputSplit = inputLine.split("\t");
                String geneID = inputSplit[0];

                if (!affinityMap.containsKey(geneID))
                {
                    affinityMap.put(geneID, new HashMap<>());
                }

                for (int i = 1; i < headerSplit.length; i++)
                {
                    String tf = headerSplit[i];

                    if (!affinityMap.get(geneID).containsKey(tf))
                    {
                        affinityMap.get(geneID).put(tf, 0.0);
                    }

                    affinityMap.get(geneID)
                            .put(tf, affinityMap.get(geneID).get(tf) + Double.parseDouble(inputSplit[i]));
                }
            }
        } catch (IOException e)
        {
            e.printStackTrace();
        }
    }

    private void computeMeanAffinities(Map<String, Map<String, Double>> affinityMap, int count)
    {
        for (String geneID : affinityMap.keySet())
        {
            for (String tf : affinityMap.get(geneID).keySet())
            {
                affinityMap.get(geneID).put(tf, affinityMap.get(geneID).get(tf) / count);
            }
        }
    }

    private void calculate(File d_input1, File d_input2, File f_output1, File f_output2, File f_ratios)
    {
        System.gc();
        String suffix = getFileSuffix();
        FileFilter filter = getSuffixFilter(suffix);
        Map<String, Map<String, Double>> affinities1 = new HashMap<>(), affinities2 = new HashMap<>();

        int count1 = Objects.requireNonNull(d_input1.listFiles(filter)).length;
        int count2 = Objects.requireNonNull(d_input2.listFiles(filter)).length;

        String header = null;

        for (File f_input : Objects.requireNonNull(d_input1.listFiles(filter)))
        {
            if (header == null)
            {
                try (BufferedReader reader = new BufferedReader(new FileReader(f_input)))
                {
                    header = reader.readLine();
                } catch (IOException e)
                {
                    logger.error(e.getMessage());
                }
                assert header != null;
            }
            mergeAffinityFile(f_input, affinities1, header);
        }
        assert header != null;

        for (File f_input : Objects.requireNonNull(d_input2.listFiles(filter)))
        {
            mergeAffinityFile(f_input, affinities2, header);
        }

        Set<String> allGenes = new HashSet<>()
        {{
            addAll(affinities1.keySet());
            addAll(affinities2.keySet());
        }};

        List<String> tfs = Arrays.asList(header.split("\t")).subList(1, header.split("\t").length);

        for (String gene : allGenes)
        {
            if (!affinities1.containsKey(gene))
            {
                affinities1.put(gene, new HashMap<>());

                for (String tf : tfs)
                {
                    affinities1.get(gene).put(tf, 0.0);
                }
            }
            if (!affinities2.containsKey(gene))
            {
                affinities2.put(gene, new HashMap<>());

                for (String tf : tfs)
                {
                    affinities2.get(gene).put(tf, 0.0);
                }
            }
        }

        logger.debug("Calculating mean affinities for group 1");
        computeMeanAffinities(affinities1, count1);
        logger.debug("Calculating mean affinities for group 2");
        computeMeanAffinities(affinities2, count2);

        Map<String, Map<String, Double>> ratios = new HashMap<>();

        logger.debug("Calculating ratios");

        for (String gene : allGenes)
        {
            ratios.put(gene, new HashMap<>());
            for (String tf : tfs)
            {
                ratios.get(gene).put(tf, (affinities1.get(gene).get(tf) + 1) / (affinities2.get(gene).get(tf) + 1));
            }
        }

        affinityMapToFile(affinities1, header, tfs, f_output1);
        affinityMapToFile(affinities2, header, tfs, f_output2);
        affinityMapToFile(ratios, header, tfs, f_ratios);
    }

    private synchronized void affinityMapToFile(Map<String, Map<String, Double>> map, String header, List<String> tfs,
                                                File target)
    {
        makeSureFileExists(target, logger);

        logger.debug("Writing affinity file: " + target.getAbsolutePath());

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(target)))
        {
            writer.write(header);
            writer.newLine();
            for (Map.Entry<String, Map<String, Double>> gene : map.entrySet())
            {
                writer.write(gene.getKey());
                for (String tf : tfs)
                {
                    writer.write("\t");
                    writer.write(String.valueOf(gene.getValue().get(tf)));
                }
                writer.newLine();
            }
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }
    }

    @Override protected int getThreadNumber()
    {
        return 2;
    }
}
