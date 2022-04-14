package lib.DistributionAnalysis;

import lib.DistributionAnalysis.Classes.TargetGeneDcg;
import lib.DistributionAnalysis.Classes.TfTargetGeneDcg;
import lib.ExecutableStep;
import tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.FileFilters.Filters;

import java.io.*;
import java.util.*;

import static util.FileManagement.extend;
import static util.FileManagement.makeSureFileExists;

public class getTopKTargetGenes extends ExecutableStep
{
    private final AbstractConfig<File> f_input_dcg = TFPRIO.configs.distributionAnalysis.fileStructure.f_dcg_stats;
    private final AbstractConfig<File> d_input_affinityValues =
            TFPRIO.configs.tepic.fileStructure.d_postprocessing_input;

    private final AbstractConfig<File> d_output = TFPRIO.configs.distributionAnalysis.fileStructure.d_dcg_targetGenes;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(Arrays.asList(f_input_dcg, d_input_affinityValues));
    }

    @Override protected Set<AbstractConfig<File>> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(d_output));
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>();
    }

    @Override protected void execute()
    {
        ArrayList<String> tfSymbols = new ArrayList<>();

        try (BufferedReader br_dcg = new BufferedReader(new FileReader(f_input_dcg.get())))
        {
            String line_dcg;
            br_dcg.readLine();
            while ((line_dcg = br_dcg.readLine()) != null)
            {
                String[] split = line_dcg.split("\t");
                tfSymbols.add(split[1]);
            }
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        for (File d_group : Objects.requireNonNull(d_input_affinityValues.get().listFiles(Filters.directoryFilter)))
        {
            for (File d_hm : Objects.requireNonNull(d_group.listFiles(Filters.directoryFilter)))
            {
                HashMap<String, TfTargetGeneDcg> tfSymbol_targetGeneDcg = new HashMap<>();
                for (String tfSymbol : tfSymbols)
                {
                    tfSymbol_targetGeneDcg.put(tfSymbol, new TfTargetGeneDcg(tfSymbol));
                }

                for (File f_sample : Objects.requireNonNull(d_hm.listFiles(Filters.fileFilter)))
                {
                    HashMap<Integer, String> index_tfSymbol = new HashMap<>();

                    try (BufferedReader reader = new BufferedReader(new FileReader(f_sample)))
                    {
                        String line = reader.readLine();
                        String[] split_header = line.split("\t");
                        for (int i = 1; i < split_header.length; i++)
                        {
                            String name_tf = split_header[i].split("_")[0];
                            name_tf = name_tf.replaceAll(":", "\\.");

                            if (tfSymbols.contains(name_tf.toUpperCase()))
                            {
                                index_tfSymbol.put(i, name_tf.toUpperCase());
                            }
                        }

                        while ((line = reader.readLine()) != null)
                        {
                            String[] split = line.split("\t");
                            String geneID = split[0];
                            String geneSymbol;

                            try
                            {
                                geneSymbol = TFPRIO.mapSymbolAndEnsg.ensgToSymbol(geneID);

                                for (int i = 1; i < split.length; i++)
                                {
                                    if (index_tfSymbol.containsKey(i))
                                    {
                                        String tfSymbol = index_tfSymbol.get(i);

                                        TfTargetGeneDcg tf_dcg = tfSymbol_targetGeneDcg.get(tfSymbol);

                                        tf_dcg.addAffinityValue(geneSymbol, Double.parseDouble(split[i]));
                                    }
                                }
                            } catch (NoSuchFieldException ignore)
                            {
                            }
                        }
                    } catch (IOException e)
                    {
                        logger.error(e.getMessage());
                    }
                }

                for (String tfSymbol : tfSymbol_targetGeneDcg.keySet())
                {
                    TfTargetGeneDcg tf_dcg = tfSymbol_targetGeneDcg.get(tfSymbol);

                    if (tf_dcg.isEmpty())
                    {
                        continue;
                    }

                    File targetFile = extend(d_output.get(), d_group.getName(), d_hm.getName(), tfSymbol + ".tsv");
                    makeSureFileExists(targetFile, logger);
                    try (BufferedWriter writer = new BufferedWriter(new FileWriter(targetFile)))
                    {
                        writer.write("TARGET_GENE\tAFFINITY_SCORE");
                        writer.newLine();
                        for (TargetGeneDcg tg_dcg : tf_dcg.getOrderedTargetGeneList())
                        {
                            writer.write(tg_dcg.toString());
                            writer.newLine();
                        }
                    } catch (IOException e)
                    {
                        logger.error(e.getMessage());
                    }
                }
            }
        }
    }
}
