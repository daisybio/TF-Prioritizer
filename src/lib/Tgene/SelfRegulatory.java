package lib.Tgene;

import lib.BinaryTree.ChromosomeRegionTrees;
import lib.ExecutableStep;
import lib.GeneRegion;
import lib.Region;
import tfprio.TFPRIO;
import util.Configs.Config;

import java.io.*;
import java.util.*;

import static util.FileManagement.extend;
import static util.FileManagement.makeSureFileExists;

public class SelfRegulatory extends ExecutableStep
{
    private final Config<File> d_input_tgene = TFPRIO.configs.tgene.fileStructure.d_groups;
    private final Config<File> d_input_tepic = TFPRIO.configs.tepic.fileStructure.d_postprocessing_output;
    private final Config<File> f_inputRegions = TFPRIO.configs.tgene.fileStructure.f_preprocessing_regions;

    private final Config<File> d_output = TFPRIO.configs.tgene.fileStructure.d_integrate;

    private final Config<String> s_tgene_groups = TFPRIO.configs.tgene.fileStructure.s_groups_mergedGroups;
    private final Config<String> s_tepic_ratios = TFPRIO.configs.tepic.fileStructure.s_postprocessing_output_ratiosDir;
    private final Config<String> consensusCalc = TFPRIO.configs.tgene.consensusCalc;
    private final Config<Double> consensus = TFPRIO.configs.tgene.consensus;

    @Override protected Set<Config<File>> getRequiredFileStructure()
    {
        return new HashSet<>(Arrays.asList(d_input_tgene, d_input_tepic, f_inputRegions));
    }

    @Override protected Set<Config<File>> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(d_output));
    }

    @Override protected Set<Config<?>> getRequiredConfigs()
    {
        return new HashSet<>(Arrays.asList(s_tgene_groups, s_tepic_ratios, consensusCalc, consensus));
    }

    @Override protected void execute()
    {
        //integrate data between TEPIC and TGENE (modifying quotients)
        //based on structure build necessary files

        for (Map.Entry<String, Set<String>> groupPairingEntry : TFPRIO.groupCombinationsToHms.entrySet())
        {
            String pairing = groupPairingEntry.getKey();

            for (String hm : groupPairingEntry.getValue())
            {
                executorService.submit(() ->
                {
                    File targetFile = extend(d_output.get(), pairing, hm, pairing + ".txt");
                    makeSureFileExists(targetFile, logger);

                    try (BufferedWriter writer = new BufferedWriter(new FileWriter(targetFile)))
                    {
                        File input_data_TGENE = extend(d_input_tgene.get(), pairing, hm, s_tgene_groups.get());
                        File input_data_TEPIC =
                                extend(d_input_tepic.get(), pairing, hm, s_tepic_ratios.get(), pairing + ".txt");

                        List<Region> preprocessedRegions = getGeneRegions(f_inputRegions.get());
                        List<Region> tgeneRegions = getGeneRegions(input_data_TGENE);

                        ChromosomeRegionTrees trees = new ChromosomeRegionTrees();
                        trees.addAllOptimized(preprocessedRegions);

                        Map<String, Set<String>> geneID_tf = new HashMap<>();

                        for (Region region : tgeneRegions)
                        {
                            GeneRegion geneRegion = (GeneRegion) region;
                            GeneRegion match = (GeneRegion) trees.getMatchingChild(region);
                            if (match != null)
                            {
                                if (!geneID_tf.containsKey(match.getId()))
                                {
                                    geneID_tf.put(match.getId(), new HashSet<>());
                                }
                                geneID_tf.get(match.getId()).add(geneRegion.getId());
                            }
                        }

                        try (BufferedReader br_tepic = new BufferedReader(new FileReader(input_data_TEPIC)))
                        {
                            String header = br_tepic.readLine();

                            writer.write(header);
                            writer.newLine();

                            HashMap<String, Integer> tfIndex = new HashMap<>();

                            ArrayList<String> tfs = new ArrayList<>(Arrays.asList(header.split("\t")));

                            for (int i = 0; i < tfs.size(); i++)
                            {
                                String key = tfs.get(i).split("_")[0].toUpperCase();
                                tfs.set(i, key);
                                String[] split_doubles = key.split("::");
                                for (String split_double : split_doubles)
                                {
                                    tfIndex.put(split_double, i);
                                }
                            }
                            String line_tepic;

                            while ((line_tepic = br_tepic.readLine()) != null)
                            {
                                String[] split = line_tepic.split("\t");

                                StringBuilder sb = new StringBuilder();
                                sb.append(split[0]);

                                if (geneID_tf.containsKey(split[0]))
                                {
                                    Set<String> tgene_pref_tfs = geneID_tf.get(split[0]);
                                    Set<Integer> pref_positions = new HashSet<>();

                                    for (String tf : tgene_pref_tfs)
                                    {
                                        pref_positions.add(tfIndex.get(tf));
                                    }

                                    for (int i = 1; i < split.length; i++)
                                    {
                                        double score = Double.parseDouble(split[i]);
                                        if (pref_positions.contains(i))
                                        {
                                            if (consensusCalc.get().equals("INCREASE_TGENE_TFS"))
                                            {
                                                score = score + (score * (1 - consensus.get()));

                                            }
                                        } else
                                        {
                                            if (consensusCalc.get().equals("DECREASE_NOT_TGENE_TFs"))
                                            {
                                                score *= (1 - consensus.get());
                                            }
                                        }
                                        sb.append("\t");
                                        sb.append(score);
                                    }
                                } else
                                {
                                    for (int i = 1; i < split.length; i++)
                                    {
                                        double score = Double.parseDouble(split[i]);
                                        if (consensusCalc.get().equals("DECREASE_NOT_TGENE_TFs"))
                                        {
                                            score *= (1 - consensus.get());
                                        }
                                        sb.append("\t");
                                        sb.append(score);
                                    }
                                }

                                writer.write(sb.toString());
                                writer.newLine();
                            }
                        }
                    } catch (IOException e)
                    {
                        logger.error(e.getMessage());
                    }
                });
            }
        }
    }

    private List<Region> getGeneRegions(File input)
    {
        List<Region> regions = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(input)))
        {
            String inputLine;
            boolean headingPossible = true;

            while ((inputLine = reader.readLine()) != null)
            {
                String[] split = inputLine.split("\t");

                if (headingPossible)
                {
                    headingPossible = false;

                    if (split[0].equals("CHR"))
                    {
                        continue;
                    }
                }

                String chromosome = split[0];
                int start = Integer.parseInt(split[1]);
                int end = Integer.parseInt(split[2]);
                String geneID = split[3];

                regions.add(new GeneRegion(chromosome, start, end, geneID));
            }
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        return regions;
    }
}
