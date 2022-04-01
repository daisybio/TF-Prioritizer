package lib.Plots;

import lib.ExecutableStep;
import lib.GeneAffinity;
import tfprio.TFPRIO;
import util.Configs.Config;
import util.FileFilters.Filters;

import java.io.*;
import java.util.*;

import static util.FileManagement.extend;
import static util.FileManagement.makeSureFileExists;

public class TopKTargetGenes extends ExecutableStep
{
    private final Config<File> d_input = TFPRIO.configs.plots.fileStructure.d_data;
    private final Config<File> d_inputTargetGenes = TFPRIO.configs.tepic.fileStructure.d_postprocessing_output;
    private final Config<File> d_inputTepicRaw = TFPRIO.configs.tepic.fileStructure.d_outputRaw;

    private final Config<File> d_output = TFPRIO.configs.plots.fileStructure.d_targetGenes;

    private final Config<String> s_plotData_different = TFPRIO.configs.plots.fileStructure.s_data_hmLevelDifferent;
    private final Config<String> s_plotData_same = TFPRIO.configs.plots.fileStructure.s_data_hmLevelSame;
    private final Config<Integer> tpmCutoff = TFPRIO.configs.tepic.tpmCutoff;
    private final Config<List<Double>> thresholds = TFPRIO.configs.plots.thresholds;
    private final Config<Integer> topKGenes = TFPRIO.configs.plots.topKGenes;

    @Override protected Set<Config<File>> getRequiredFileStructure()
    {
        return new HashSet<>()
        {{
            add(d_input);
            if (tpmCutoff.get() > 0)
            {
                add(d_inputTargetGenes);
            } else
            {
                add(d_inputTepicRaw);
            }
        }};
    }

    @Override protected Set<Config<File>> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(d_output));
    }

    @Override protected Set<Config<?>> getRequiredConfigs()
    {
        return new HashSet<>(Arrays.asList(s_plotData_different, s_plotData_same, tpmCutoff, thresholds, topKGenes));
    }

    @Override protected void execute()
    {
        for (String hm : TFPRIO.existingHms)
        {
            for (double threshold : thresholds.get())
            {
                for (String suffix : Arrays.asList(s_plotData_different.get(), s_plotData_same.get()))
                {
                    File f_input = extend(d_input.get(), hm, String.valueOf(threshold), suffix);
                    File d_out = extend(d_output.get(), hm, String.valueOf(threshold),
                            suffix.substring(0, suffix.lastIndexOf(".")));

                    try (BufferedReader reader = new BufferedReader(new FileReader(f_input)))
                    {
                        String line_tf = reader.readLine();
                        String[] split_header = line_tf.split(",");
                        List<String> pairingsInHeader = new ArrayList<>();

                        for (int i = 1; i < split_header.length; i++)
                        {
                            String pairing = split_header[i].substring(split_header[i].indexOf(":") + 1);
                            pairingsInHeader.add(pairing);
                        }

                        while ((line_tf = reader.readLine()) != null)
                        {
                            String[] split = line_tf.split(",");
                            String tfSymbol = split[0];

                            ArrayList<String> identifiedPairings = new ArrayList<>();

                            for (int i = 1; i < split.length; i++)
                            {
                                if (!split[i].isEmpty())
                                {
                                    identifiedPairings.add(pairingsInHeader.get(i - 1));
                                }
                            }

                            for (String pairing : identifiedPairings)
                            {
                                String[] splitPairing = pairing.split("_");
                                String group1 = splitPairing[0];
                                String group2 = splitPairing[1];

                                executorService.submit(() -> writeTargetGenes(pairing, group1, hm, d_out, tfSymbol));
                                executorService.submit(() -> writeTargetGenes(pairing, group2, hm, d_out, tfSymbol));
                            }
                        }
                    } catch (IOException e)
                    {
                        logger.error(e.getMessage());
                    }
                }
            }
        }
    }

    private void writeTargetGenes(String pairing, String group, String hm, File d_out, String tfSymbol)
    {
        File f_output = extend(d_out, group);

        File parent;
        String suffix;

        List<File> inputFiles;
        if (tpmCutoff.get() > 0)
        {
            suffix = "_Gene_View_Filtered_TPM.txt";
            parent = extend(d_inputTargetGenes.get(), pairing, hm, group);
            inputFiles =
                    new ArrayList<>(List.of(Objects.requireNonNull(parent.listFiles(Filters.getSuffixFilter(suffix)))));
        } else
        {
            suffix = "_Gene_View_Filtered.txt";
            parent = extend(d_inputTepicRaw.get(), group, hm);
            inputFiles = new ArrayList<>();

            for (File d_sub : Objects.requireNonNull(parent.listFiles(Filters.directoryFilter)))
            {
                inputFiles.addAll(List.of(Objects.requireNonNull(d_sub.listFiles(Filters.getSuffixFilter(suffix)))));
            }
        }

        for (File f_input : inputFiles)
        {
            ArrayList<GeneAffinity> all_affinities = new ArrayList<>();
            try (BufferedReader reader = new BufferedReader(new FileReader(f_input)))
            {
                String line = reader.readLine();

                Integer tfColumn = null;

                String[] header = line.split("\t");
                for (int i = 0; i < header.length; i++)
                {
                    if (header[i].toUpperCase().contains(tfSymbol.toUpperCase()))
                    {
                        tfColumn = i;
                    }
                }

                if (tfColumn == null)
                {
                    return;
                }

                while ((line = reader.readLine()) != null)
                {
                    String[] split = line.split("\t");

                    String geneID = split[0];
                    double affinityValue = Double.parseDouble(split[tfColumn]);

                    GeneAffinity affinity = new GeneAffinity(geneID, affinityValue);
                    all_affinities.add(affinity);
                }
            } catch (IOException e)
            {
                logger.error(e.getMessage());
            }

            Collections.sort(all_affinities);

            File targetFile = extend(f_output, tfSymbol + ".csv");
            makeSureFileExists(targetFile, logger);

            try (BufferedWriter writer = new BufferedWriter(new FileWriter(targetFile)))
            {
                writer.write("ENSG\tSYMBOL\tAFFINITY");
                writer.newLine();

                for (int i = 0; i < topKGenes.get(); i++)
                {
                    if (all_affinities.size() > i)
                    {
                        writer.write(all_affinities.get(i).toString());
                        writer.newLine();
                    }
                }
            } catch (IOException e)
            {
                logger.error(e.getMessage());
            }
        }
    }
}