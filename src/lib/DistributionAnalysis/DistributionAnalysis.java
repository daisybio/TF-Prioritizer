package lib.DistributionAnalysis;

import lib.ExecutableStep;
import tfprio.tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.FileFilters.Filters;

import java.io.*;
import java.util.*;

import static util.FileManagement.*;

public class DistributionAnalysis extends ExecutableStep
{
    private final AbstractConfig<File> f_inputTfs = TFPRIO.configs.distributionAnalysis.fileStructure.f_analyzedTfs;
    private final AbstractConfig<File> d_inputGeneCounts =
            TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_meanCounts;
    private final AbstractConfig<File> d_inputDifferentialExpression = TFPRIO.configs.deSeq2.fileStructure.d_output;
    private final AbstractConfig<File> d_inputTgene = TFPRIO.configs.tgene.fileStructure.d_filteredTargetGenes;
    private final AbstractConfig<File> d_inputDynamite = TFPRIO.configs.dynamite.fileStructure.d_output;
    private final AbstractConfig<File> d_inputTepicProcessed =
            TFPRIO.configs.tepic.fileStructure.d_postprocessing_output;
    private final AbstractConfig<File> d_inputTepicRaw = TFPRIO.configs.tepic.fileStructure.d_outputRaw;

    private final GeneratedFileStructure d_outputRaw =
            TFPRIO.configs.distributionAnalysis.fileStructure.d_tfTgScores_raw;
    private final GeneratedFileStructure d_outputBackgroundDistributionAll =
            TFPRIO.configs.distributionAnalysis.fileStructure.d_tfTgScores_backgroundDistribution_all;
    private final GeneratedFileStructure d_outputBackgroundDistributionHm =
            TFPRIO.configs.distributionAnalysis.fileStructure.d_tfTgScores_backgroundDistribution_hm;
    private final GeneratedFileStructure d_outputTFScoresAll =
            TFPRIO.configs.distributionAnalysis.fileStructure.d_tfTgScores_tfDistribution_all;
    private final GeneratedFileStructure d_outputTFScoresHm =
            TFPRIO.configs.distributionAnalysis.fileStructure.d_tfTgScores_tfDistribution_hm;

    private final AbstractConfig<Boolean> performAllAnalysis = TFPRIO.configs.distributionAnalysis.performAllAnalysis;
    private final AbstractConfig<String> s_dynamiteInput = TFPRIO.configs.dynamite.fileStructure.s_output_toBePlotted;
    private final AbstractConfig<Boolean> randomizeTfGeneMatrix = TFPRIO.configs.tepic.randomizeTfGeneMatrix;
    private final AbstractConfig<String> scoreType = TFPRIO.configs.plots.distributionAnalysisScoreType;
    private final AbstractConfig<String> s_output =
            TFPRIO.configs.distributionAnalysis.fileStructure.s_tfTgScores_backgroundDistribution_csv;

    // Optional configs
    private final AbstractConfig<File> f_tgeneExecutable = TFPRIO.configs.tgene.pathToExecutable;
    private final AbstractConfig<Double> tpmCutoff = TFPRIO.configs.tepic.tpmCutoff;


    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>()
        {{
            add(f_inputTfs);
            add(d_inputGeneCounts);
            add(d_inputDifferentialExpression);
            add(d_inputTgene);
            add(d_inputDynamite);
            if (tpmCutoff.isSet())
            {
                add(d_inputTepicProcessed);
            } else
            {
                add(d_inputTepicRaw);
            }
        }};
    }

    @Override public Set<GeneratedFileStructure> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(d_outputBackgroundDistributionHm, d_outputTFScoresHm, d_outputRaw))
        {{
            if (performAllAnalysis.get())
            {
                add(d_outputBackgroundDistributionAll);
                add(d_outputTFScoresAll);
            } else
            {
                for (GeneratedFileStructure fileStructure : Arrays.asList(d_outputBackgroundDistributionAll,
                        d_outputTFScoresAll))
                {
                    fileStructure.deleteAndSetNoGenerationReason(performAllAnalysis.getName() + " is set to false");
                }
            }
        }};
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>(
                Arrays.asList(performAllAnalysis, s_dynamiteInput, randomizeTfGeneMatrix, scoreType, s_output));
    }

    @Override protected Set<AbstractConfig<?>> getOptionalConfigs()
    {
        return new HashSet<>(List.of(f_tgeneExecutable, tpmCutoff));
    }

    @Override protected void execute()
    {
        HashMap<String, HashSet<String>> distinct_hms_tfs = new HashMap<>();
        HashMap<String, HashSet<String>> distinct_tfs_hms = new HashMap<>();

        HashMap<String, HashMap<String, Double>> group_geneID_count = new HashMap<>();
        HashMap<String, HashMap<String, Double>> groupPairing_geneID_differentialExpression = new HashMap<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(f_inputTfs.get())))
        {
            String inputLine;
            reader.readLine();
            while ((inputLine = reader.readLine()) != null)
            {
                String[] split = inputLine.split("\t");
                String tfSymbol = split[0];
                String hm = split[1];

                if (!distinct_hms_tfs.containsKey(hm))
                {
                    distinct_hms_tfs.put(hm, new HashSet<>());
                }
                if (!distinct_tfs_hms.containsKey(tfSymbol))
                {
                    distinct_tfs_hms.put(tfSymbol, new HashSet<>());
                }

                distinct_hms_tfs.get(hm).add(tfSymbol);
                distinct_tfs_hms.get(tfSymbol).add(hm);
            }
        } catch (IOException e)
        {
            logger.error("A " + e.getMessage());
        }

        Random random = new Random();

        for (String group : TFPRIO.groupsToHms.keySet())
        {
            File f_input = extend(d_inputGeneCounts.get(), group + ".tsv");

            HashMap<String, Double> ensg_geneCount = new HashMap<>();

            try (BufferedReader br_gene_counts = new BufferedReader(new FileReader(f_input)))
            {
                String inputLine;
                br_gene_counts.readLine();
                while ((inputLine = br_gene_counts.readLine()) != null)
                {
                    String[] split = inputLine.split("\t");
                    ensg_geneCount.put(split[1].toUpperCase(), Double.parseDouble(split[2]));
                }
            } catch (IOException e)
            {
                logger.error("B " + e.getMessage());
            }
            group_geneID_count.put(group, ensg_geneCount);
        }


        for (String pairing : TFPRIO.groupCombinationsToHms.keySet())
        {
            File f_input = extend(d_inputDifferentialExpression.get(), pairing + ".tsv");

            HashMap<String, Double> current_diff_expr = new HashMap<>();

            try (BufferedReader reader = new BufferedReader(new FileReader(f_input)))
            {
                String inputLine;
                reader.readLine();
                while ((inputLine = reader.readLine()) != null)
                {
                    String[] split = inputLine.split("\t");
                    current_diff_expr.put(split[0], Double.parseDouble(split[1]));
                }
            } catch (IOException e)
            {
                logger.error("C " + e.getMessage());
            }
            groupPairing_geneID_differentialExpression.put(pairing, current_diff_expr);
        }

        for (String tfSymbol : distinct_tfs_hms.keySet())
        {
            Set<String> geneIDs;
            try
            {
                geneIDs = TFPRIO.mapSymbolAndEnsg.symbolToEnsg(tfSymbol);
            } catch (NoSuchFieldException e)
            {
                logger.warn(e.getMessage());
                continue;
            }
            HashMap<String, Double> groupPairing_differentialExpression = new HashMap<>();

            for (String groupPairing : groupPairing_geneID_differentialExpression.keySet())
            {
                double sum = 0;
                int count = 0;

                HashMap<String, Double> geneID_differentialExpression =
                        groupPairing_geneID_differentialExpression.get(groupPairing);

                for (String geneID : geneIDs)
                {
                    if (geneID_differentialExpression.containsKey(geneID))
                    {
                        sum += geneID_differentialExpression.get(geneID);
                        count++;
                    }
                }

                groupPairing_differentialExpression.put(groupPairing, sum / count);
            }

            for (String groupPairing : groupPairing_differentialExpression.keySet())
            {
                String[] pairingSplit = groupPairing.split("_");
                String group1 = pairingSplit[0];
                String group2 = pairingSplit[1];

                for (String hm : distinct_tfs_hms.get(tfSymbol))
                {
                    if (!TFPRIO.groupCombinationsToHms.get(groupPairing).contains(hm))
                    {
                        continue;
                    }

                    executorService.submit(() ->
                    {
                        String tfSymbolExtension = "";
                        HashSet<String> available_ensgs = new HashSet<>();

                        if (f_tgeneExecutable.isSet())
                        {
                            File f_input = extend(d_inputTgene.get(), groupPairing, hm, groupPairing + ".txt");

                            try (BufferedReader reader = new BufferedReader(new FileReader(f_input)))
                            {
                                String inputLine;
                                reader.readLine();
                                while ((inputLine = reader.readLine()) != null)
                                {
                                    available_ensgs.add(inputLine.substring(0, inputLine.indexOf("\t")));
                                }
                            } catch (IOException e)
                            {
                                logger.error("D " + e.getMessage());
                            }
                        }

                        Double tf_regression_coefficient = null;

                        File f_inputDynamite = extend(d_inputDynamite.get(), groupPairing, hm, s_dynamiteInput.get());

                        try (BufferedReader reader = new BufferedReader(new FileReader(f_inputDynamite)))
                        {
                            String inputLine;
                            reader.readLine();
                            while ((inputLine = reader.readLine()) != null)
                            {
                                String[] split = inputLine.split("\t");

                                String[] tfSplit = split[0].split("_");

                                if (tfSplit[0].equalsIgnoreCase(tfSymbol))
                                {
                                    if (tfSplit.length > 1)
                                    {
                                        tfSymbolExtension = tfSplit[1];
                                    }

                                    if (!split[1].equals("0"))
                                    {
                                        tf_regression_coefficient = Double.parseDouble(split[1]);
                                    }
                                    break;
                                }
                            }
                        } catch (IOException e)
                        {
                            logger.error("E " + e.getMessage());
                        }

                        if (tf_regression_coefficient != null)
                        {
                            List<File> inputFiles1 = getTepicInputFiles(groupPairing, hm, group1);
                            List<File> inputFiles2 = getTepicInputFiles(groupPairing, hm, group2);

                            String tfSymbolManipulated = tfSymbol.replace(".", ":");
                            String tfSymbolExtended = tfSymbol.replace(".", ":") + "_" + tfSymbolExtension;

                            List<String> acceptedTfNames =
                                    Arrays.asList(tfSymbolManipulated, tfSymbolExtended, tfSymbol);

                            Map<String, Double> geneID_score_1 = averageGeneScores(inputFiles1, acceptedTfNames);
                            Map<String, Double> geneID_score_2 = averageGeneScores(inputFiles2, acceptedTfNames);

                            File targetFile = extend(d_outputRaw.get(), hm, tfSymbol, groupPairing + ".tsv");
                            makeSureFileExists(targetFile, logger);

                            try (BufferedWriter writer = new BufferedWriter(new FileWriter(targetFile)))
                            {
                                for (String geneID : geneID_score_1.keySet())
                                {
                                    if (!geneID_score_2.containsKey(geneID) ||
                                            (f_tgeneExecutable.isSet() && !available_ensgs.contains(geneID)) ||
                                            !group_geneID_count.get(group1).containsKey(geneID) ||
                                            !group_geneID_count.get(group2).containsKey(geneID) ||
                                            !groupPairing_differentialExpression.containsKey(groupPairing))
                                    {
                                        continue;
                                    }

                                    double geneScore1 = geneID_score_1.get(geneID);
                                    double geneScore2 = geneID_score_2.get(geneID);

                                    double geneCount1 = group_geneID_count.get(group1).get(geneID);
                                    double geneCount2 = group_geneID_count.get(group2).get(geneID);

                                    double diffExpression;

                                    if (randomizeTfGeneMatrix.get())
                                    {
                                        List<String> differentialExpressionList =
                                                new ArrayList<>(groupPairing_differentialExpression.keySet());

                                        int randomIndex = random.nextInt(differentialExpressionList.size());
                                        String randomPairing = differentialExpressionList.get(randomIndex);
                                        diffExpression = groupPairing_differentialExpression.get(randomPairing);
                                    } else
                                    {
                                        diffExpression = groupPairing_differentialExpression.get(groupPairing);
                                    }

                                    try
                                    {
                                        double tgScore1;
                                        double tgScore2;

                                        if (scoreType.get().equals("EXCL_GENE_COUNTS"))
                                        {
                                            tgScore1 = diffExpression * geneScore1;
                                            tgScore2 = diffExpression * geneScore2;
                                        } else if (scoreType.get().equals("GENE_COUNTS"))
                                        {
                                            tgScore1 = geneCount1 * diffExpression * geneScore1;
                                            tgScore2 = geneCount2 * diffExpression * geneScore2;
                                        } else {
                                            tgScore1 = 0;
                                            tgScore2 = 0;
                                        }

                                        tgScore1 = Math.abs(tgScore1);
                                        tgScore2 = Math.abs(tgScore2);

                                        double tgScoreCombined = tgScore1 + tgScore2;

                                        double tfTgScore = Math.abs(tgScoreCombined * tf_regression_coefficient);
                                        if (!Double.isNaN(tfTgScore))
                                        {
                                            String line =
                                                    geneID + "\t" + tfTgScore + "\t" + hm + "\t" + groupPairing + "\t" +
                                                            tfSymbol + "\t" + tf_regression_coefficient;

                                            writer.write(line);
                                            writer.newLine();
                                        }
                                    } catch (NullPointerException | AssertionError e)
                                    {
                                        logger.warn(groupPairing + ", " + hm + ", " + tfSymbol + ": " + e.getMessage());
                                    }
                                }
                            } catch (IOException e)
                            {
                                logger.error(e.getMessage());
                            }
                        }
                    });
                }
            }
        }
        finishAllQueuedThreads();
        logger.info("Start concatenating files");
        concatFiles(distinct_tfs_hms.keySet());
    }

    private void concatFiles(Iterable<String> tfSymbols)
    {
        StringBuilder sb_all = new StringBuilder();

        Map<String, StringBuilder> map_sbAllHms = new HashMap<>();

        try
        {
            for (String hm : TFPRIO.existingHms)
            {
                StringBuilder sb_allTfs = new StringBuilder();

                for (String tfsymbol : tfSymbols)
                {
                    StringBuilder sb_tf = new StringBuilder();

                    boolean foundData = false;

                    for (String pairing : TFPRIO.groupCombinationsToHms.keySet())
                    {
                        File sourceFile = extend(d_outputRaw.get(), hm, tfsymbol, pairing + ".tsv");

                        if (sourceFile.exists())
                        {
                            foundData = true;
                            String content = readFile(sourceFile);

                            sb_tf.append(content);
                            sb_allTfs.append(content);

                            if (performAllAnalysis.get())
                            {
                                sb_all.append(content);

                                if (!map_sbAllHms.containsKey(tfsymbol))
                                {
                                    map_sbAllHms.put(tfsymbol, new StringBuilder());
                                }
                                map_sbAllHms.get(tfsymbol).append(content);
                            }
                        }
                    }

                    if (foundData)
                    {
                        File f_tf = extend(d_outputTFScoresHm.get(), hm, tfsymbol + ".tsv");
                        if (sb_tf.length() > 0)
                        {
                            writeFile(f_tf, getHeader(hm, tfsymbol) + sb_tf);
                        }
                    }
                }

                File f_allTfs = extend(d_outputBackgroundDistributionHm.get(), hm, s_output.get());
                writeFile(f_allTfs, getHeader(hm, "ALL") + sb_allTfs);
            }

            if (performAllAnalysis.get())
            {
                File f_all = extend(d_outputBackgroundDistributionAll.get(), s_output.get());
                if (sb_all.length() > 0)
                {
                    writeFile(f_all, getHeader("ALL", "ALL") + sb_all);
                }

                for (Map.Entry<String, StringBuilder> allHmsEntry : map_sbAllHms.entrySet())
                {
                    File f_allHms = extend(d_outputTFScoresAll.get(), allHmsEntry.getKey() + ".tsv");
                    if (allHmsEntry.getValue().length() > 0)
                    {
                        writeFile(f_allHms, getHeader("ALL", allHmsEntry.getKey()) + allHmsEntry.getValue());
                    }
                }
            }
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }
    }

    private List<File> getTepicInputFiles(String pairing, String hm, String group)
    {
        List<File> inputFiles;
        String suffix;
        if (tpmCutoff.isSet())
        {
            suffix = "_Gene_View_Filtered_TPM.txt";

            File parent = extend(d_inputTepicProcessed.get(), pairing, hm, group);
            inputFiles =
                    new ArrayList<>(List.of(Objects.requireNonNull(parent.listFiles(Filters.getSuffixFilter(suffix)))));

        } else
        {
            suffix = "_Gene_View_Filtered.txt";

            File parent = extend(d_inputTepicRaw.get(), group, hm);
            inputFiles = new ArrayList<>();

            for (File d_sub : Objects.requireNonNull(parent.listFiles(Filters.directoryFilter)))
            {
                inputFiles.addAll(List.of(Objects.requireNonNull(d_sub.listFiles(Filters.getSuffixFilter(suffix)))));
            }
        }

        return inputFiles;
    }

    private Map<String, Double> averageGeneScores(List<File> inputFiles, List<String> tfSymbols)
    {
        Map<String, Double> scores = new HashMap<>();
        Map<String, Integer> counts = new HashMap<>();

        for (File inputFile : inputFiles)
        {
            try (BufferedReader reader = new BufferedReader(new FileReader(inputFile)))
            {
                String header = reader.readLine();

                Integer tfColumn = null;
                String[] headerSplit = header.split("\t");
                for (int i = 0; i < headerSplit.length; i++)
                {
                    if (tfSymbols.contains(headerSplit[i].toUpperCase()))
                    {
                        tfColumn = i;
                        break;
                    }
                }

                if (tfColumn == null)
                {
                    continue;
                }

                String inputLine;
                while ((inputLine = reader.readLine()) != null)
                {
                    String[] split = inputLine.split("\t");
                    String geneID = split[0];
                    double value = Double.parseDouble(split[tfColumn]);

                    if (!scores.containsKey(geneID))
                    {
                        scores.put(geneID, .0);
                    }
                    if (!counts.containsKey(geneID))
                    {
                        counts.put(geneID, 0);
                    }
                    scores.put(geneID, scores.get(geneID) + value);
                    counts.put(geneID, counts.get(geneID) + 1);
                }
            } catch (IOException e)
            {
                logger.error("G " + e.getMessage());
            }
        }

        for (String geneID : scores.keySet())
        {
            if (!counts.containsKey(geneID))
            {
                logger.warn("Unrecognized geneID: " + geneID);
                continue;
            }

            scores.put(geneID, scores.get(geneID) / counts.get(geneID));
        }

        return scores;
    }

    private String getHeader(String hm, String tf)
    {
        return "##HM\t" + hm + "\n" + "##TF\t" + tf + "\n" + "TARGET_GENE\tTF_TG_SCORE\tHM\tGROUPS\tTF\tTF_COEFF\n";
    }

    @Override protected int getMemoryEstimationMb()
    {
        return 5000;
    }
}
