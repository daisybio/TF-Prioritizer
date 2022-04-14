package lib.DistributionAnalysis;

import lib.ExecutableStep;
import tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.FileFilters.Filters;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

import static util.FileManagement.makeSureFileExists;
import static util.ScriptExecution.executeAndWait;

public class CoOccurrenceAnalysis extends ExecutableStep
{
    private final AbstractConfig<File> d_input = TFPRIO.configs.chipAtlas.fileStructure.d_peakFiles;

    private final AbstractConfig<File> f_output_concatenated =
            TFPRIO.configs.distributionAnalysis.fileStructure.f_cooccurrence_concatenatedBed;
    private final AbstractConfig<File> f_output_merged =
            TFPRIO.configs.distributionAnalysis.fileStructure.f_cooccurrence_mergedBed;
    private final AbstractConfig<File> f_output_sorted =
            TFPRIO.configs.distributionAnalysis.fileStructure.f_cooccurrence_sortedBed;
    private final AbstractConfig<File> f_output_frequencies =
            TFPRIO.configs.distributionAnalysis.fileStructure.f_cooccurrence_frequencies;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>();
    }

    @Override protected Set<AbstractConfig<File>> getCreatedFileStructure()
    {
        return new HashSet<>();
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>();
    }

    @Override protected void execute()
    {
        DecimalFormat df = new DecimalFormat("0.00");

        makeSureFileExists(f_output_concatenated.get(), logger);
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(f_output_concatenated.get())))
        {
            for (File d_tf : Objects.requireNonNull(d_input.get().listFiles(Filters.directoryFilter)))
            {
                for (File f_bedFile : Objects.requireNonNull(d_tf.listFiles(Filters.getSuffixFilter(".bed"))))
                {
                    String name = f_bedFile.getName().substring(0, f_bedFile.getName().lastIndexOf("."));

                    try (BufferedReader reader = new BufferedReader(new FileReader(f_bedFile)))
                    {
                        String inputLine;
                        while ((inputLine = reader.readLine()) != null)
                        {
                            if (inputLine.startsWith("track"))
                            {
                                continue;
                            }

                            String[] split = inputLine.split("\t");

                            writer.write(split[0]);
                            writer.write("\t");
                            writer.write(split[1]);
                            writer.write("\t");
                            writer.write(split[2]);
                            writer.write("\t");
                            writer.write(name);
                            writer.newLine();
                        }
                    }
                }
            }
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        executeAndWait("sort -k1,1 -k2,2n " + f_output_concatenated.get().getAbsolutePath() + " -o " +
                f_output_sorted.get().getAbsolutePath(), logger);

        executeAndWait(Arrays.asList("/bin/bash", "-c",
                "bedtools merge -i " + f_output_sorted.get().getAbsolutePath() + " -c 4 -o " + "collapse -delim " +
                        "'|' > " + f_output_merged.get().getAbsolutePath()), logger);

        HashMap<String, Integer> tf_total_occurrences = new HashMap<>();
        HashMap<String, HashMap<String, Integer>> tf_otherTf_cooccurrenceCount = new HashMap<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(f_output_merged.get())))
        {
            String inputLine;
            while ((inputLine = reader.readLine()) != null)
            {
                String[] split = inputLine.split("\t");

                if (split.length != 4)
                {
                    continue;
                }

                List<String> tfSymbols = Arrays.asList(split[3].split("\\|"));

                for (String tfSymbol : tfSymbols)
                {
                    if (tf_total_occurrences.containsKey(tfSymbol))
                    {
                        tf_total_occurrences.put(tfSymbol, tf_total_occurrences.get(tfSymbol) + 1);
                    } else
                    {
                        tf_total_occurrences.put(tfSymbol, 1);
                    }
                }

                for (String tfSymbol : tfSymbols)
                {
                    if (!tf_otherTf_cooccurrenceCount.containsKey(tfSymbol))
                    {
                        tf_otherTf_cooccurrenceCount.put(tfSymbol, new HashMap<>());
                    }

                    for (String otherTfSymbol : tfSymbols)
                    {
                        if (tfSymbol.equals(otherTfSymbol))
                        {
                            continue;
                        }

                        if (!tf_otherTf_cooccurrenceCount.get(tfSymbol).containsKey(otherTfSymbol))
                        {
                            tf_otherTf_cooccurrenceCount.get(tfSymbol).put(otherTfSymbol, 0);
                        }

                        tf_otherTf_cooccurrenceCount.get(tfSymbol)
                                .put(otherTfSymbol, tf_otherTf_cooccurrenceCount.get(tfSymbol).get(otherTfSymbol) + 1);
                    }
                }
            }
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        HashMap<String, HashMap<String, Double>> tf_otherTf_cooccurrencePercentage = new HashMap<>();
        for (String tfSymbol : tf_otherTf_cooccurrenceCount.keySet())
        {
            HashMap<String, Double> otherTf_cooccurrencePercentage = new HashMap<>();
            tf_otherTf_cooccurrencePercentage.put(tfSymbol, otherTf_cooccurrencePercentage);

            HashMap<String, Integer> otherTf_counts = tf_otherTf_cooccurrenceCount.get(tfSymbol);
            for (String otherTfSymbol : otherTf_counts.keySet())
            {
                int total_counts = tf_total_occurrences.get(tfSymbol) + tf_total_occurrences.get(otherTfSymbol);
                int co_occurences = otherTf_counts.get(otherTfSymbol);

                double freq = (co_occurences * 1.0) / total_counts;

                otherTf_cooccurrencePercentage.put(otherTfSymbol, freq);
            }
        }

        makeSureFileExists(f_output_frequencies.get(), logger);
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(f_output_frequencies.get())))
        {
            ArrayList<String> tfSymbols = new ArrayList<>(tf_otherTf_cooccurrencePercentage.keySet());

            for (String tfSymbol : tfSymbols)
            {
                writer.write("\t");
                writer.write(tfSymbol);
            }
            writer.newLine();

            for (String tfSymbol : tfSymbols)
            {
                writer.write(tfSymbol);
                for (String otherTfSymbol : tfSymbols)
                {
                    writer.write("\t");
                    if (tfSymbol.equals(otherTfSymbol))
                    {
                        writer.write("1.00");
                    } else
                    {
                        double freq_column_row = tf_otherTf_cooccurrencePercentage.get(otherTfSymbol).get(tfSymbol);
                        writer.write(df.format(freq_column_row));
                    }
                }
                writer.newLine();
            }
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }
    }
}
