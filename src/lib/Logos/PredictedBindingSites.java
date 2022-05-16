package lib.Logos;

import lib.ExecutableStep;
import tfprio.tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.FileFilters.Filters;

import java.io.*;
import java.util.*;

import static util.FileManagement.*;
import static util.ScriptExecution.executeAndWait;

public class PredictedBindingSites extends ExecutableStep
{
    private final AbstractConfig<File> f_input_dcgResult =
            TFPRIO.configs.distributionAnalysis.fileStructure.f_dcg_stats;
    private final AbstractConfig<File> d_input_sequences = TFPRIO.configs.tepic.fileStructure.d_outputRaw;
    private final AbstractConfig<File> f_scriptTemplate = TFPRIO.configs.scriptTemplates.f_logos_tfBindingSequence;

    private final GeneratedFileStructure d_output_data =
            TFPRIO.configs.distributionAnalysis.fileStructure.d_logos_tfBindingSequence_data;

    private final AbstractConfig<String> s_output_script =
            TFPRIO.configs.distributionAnalysis.fileStructure.s_logos_tfBindingSequence_script;
    private final AbstractConfig<String> s_fasta =
            TFPRIO.configs.distributionAnalysis.fileStructure.s_logos_tfBindingSequence_fasta;
    private final AbstractConfig<String> allName = TFPRIO.configs.distributionAnalysis.allName;
    private final AbstractConfig<String> s_trapSequences = TFPRIO.configs.tepic.fileStructure.s_outputRaw_trapSequences;
    private final AbstractConfig<Double> affinityCutoff = TFPRIO.configs.plots.trapPredictedSequenceLogosAffinityCutoff;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(Arrays.asList(f_input_dcgResult, d_input_sequences, f_scriptTemplate));
    }

    @Override public Set<GeneratedFileStructure> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(d_output_data));
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>(Arrays.asList(s_output_script, s_fasta, allName, s_trapSequences, affinityCutoff));
    }

    @Override protected void execute()
    {
        // TODO: Improve multithreading
        // Currently all the input files are read for all the tfs, there has to be a better solution

        String scriptTemplate = null;
        try
        {
            scriptTemplate = readFile(f_scriptTemplate.get());
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }
        assert scriptTemplate != null;

        ArrayList<String> tfSymbols = new ArrayList<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(f_input_dcgResult.get())))
        {
            String inputLine;
            reader.readLine();
            while ((inputLine = reader.readLine()) != null)
            {
                String[] split = inputLine.split("\t");
                tfSymbols.add(split[1]);
            }
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        int rank = 1;
        for (String tfSymbol : tfSymbols)
        {
            final int finalRank = rank;
            final String finalScriptTemplate = scriptTemplate;
            executorService.submit(() ->
            {
                HashMap<String, BufferedWriter> hm_writer = new HashMap<>();
                File f_output_allFasta =
                        extend(d_output_data.get(), finalRank + "_" + tfSymbol, allName.get() + s_fasta.get());
                makeSureFileExists(f_output_allFasta, logger);

                HashMap<String, Integer> bufferedWriters_sequencesNumbers = new HashMap<>();
                bufferedWriters_sequencesNumbers.put(allName.get(), 0);

                HashMap<String, File> hm_files = new HashMap<>();
                hm_files.put(allName.get(), f_output_allFasta);

                try
                {
                    hm_writer.put(allName.get(), new BufferedWriter(new FileWriter(f_output_allFasta)));
                    for (File d_group : Objects.requireNonNull(
                            d_input_sequences.get().listFiles(Filters.directoryFilter)))
                    {
                        for (File f_hm : Objects.requireNonNull(d_group.listFiles(Filters.directoryFilter)))
                        {
                            String hm = f_hm.getName();

                            if (!hm_writer.containsKey(hm))
                            {
                                File f_hmFasta =
                                        extend(d_output_data.get(), finalRank + "_" + tfSymbol, hm + s_fasta.get());
                                makeSureFileExists(f_hmFasta, logger);
                                hm_writer.put(hm, new BufferedWriter(new FileWriter(f_hmFasta)));
                                bufferedWriters_sequencesNumbers.put(hm, 0);
                                hm_files.put(hm, f_hmFasta);
                            }
                        }
                    }
                } catch (IOException e)
                {
                    logger.error(e.getMessage());
                }

                //run through all runs and add sequences to fasta
                for (File d_group : Objects.requireNonNull(d_input_sequences.get().listFiles(Filters.directoryFilter)))
                {
                    for (File d_hm : Objects.requireNonNull(d_group.listFiles(Filters.directoryFilter)))
                    {
                        String hm = d_hm.getName();

                        for (File d_sample : Objects.requireNonNull(d_hm.listFiles(Filters.directoryFilter)))
                        {
                            for (File f_sample : Objects.requireNonNull(
                                    d_sample.listFiles(Filters.getSuffixFilter(s_trapSequences.get()))))
                            {
                                try (BufferedReader reader = new BufferedReader(new FileReader(f_sample)))
                                {
                                    String inputLine;
                                    while ((inputLine = reader.readLine()) != null)
                                    {
                                        if (inputLine.startsWith("#") || inputLine.startsWith("TF\tAFFINITY_VALUE"))
                                        {
                                            continue;
                                        }

                                        String[] split = inputLine.split("\t");

                                        String name_tf = split[0].split("_")[0].toUpperCase().replace(":", ".");

                                        if (tfSymbol.equals(name_tf))
                                        {
                                            double affinity = Double.parseDouble(split[1]);
                                            if (affinity < affinityCutoff.get())
                                            {
                                                continue;
                                            }

                                            String sequence = split[5];

                                            {
                                                BufferedWriter writer = hm_writer.get(allName.get());
                                                int seqCount = bufferedWriters_sequencesNumbers.get(allName.get());
                                                writer.write(">" + seqCount);
                                                writer.newLine();
                                                writer.write(sequence);
                                                writer.newLine();
                                                seqCount++;
                                                bufferedWriters_sequencesNumbers.put(allName.get(), seqCount);
                                            }

                                            {
                                                BufferedWriter writer = hm_writer.get(hm);
                                                int seqCount = bufferedWriters_sequencesNumbers.get(hm);
                                                writer.write(">" + seqCount);
                                                writer.newLine();
                                                writer.write(sequence);
                                                writer.newLine();
                                                seqCount++;
                                                bufferedWriters_sequencesNumbers.put(hm, seqCount);
                                            }
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

                for (Map.Entry<String, BufferedWriter> entry : hm_writer.entrySet())
                {
                    try
                    {
                        entry.getValue().close();
                    } catch (IOException e)
                    {
                        logger.error(e.getMessage());
                    }
                }

                //create frequencies
                File f_output_data = extend(d_output_data.get(), finalRank + "_" + tfSymbol);
                File f_output_logo = extend(d_output_data.get(), finalRank + "_" + tfSymbol);

                StringBuilder sb_calls = new StringBuilder();

                for (String hm : hm_files.keySet())
                {
                    File f_logo = extend(f_output_logo, tfSymbol + "_" + hm + ".png");
                    File f_input = hm_files.get(hm);
                    File f_output = extend(f_input.getParentFile(), hm + ".motif");

                    sb_calls.append("generate('").append(f_input.getAbsolutePath()).append("', '")
                            .append(f_logo.getAbsolutePath()).append("', '").append(f_output).append("')\n");
                }

                File f_script = extend(f_output_data, s_output_script.get());
                String script = finalScriptTemplate.replace("{CALLS}", sb_calls.toString());

                try
                {
                    writeFile(f_script, script);
                } catch (IOException e)
                {
                    logger.error(e.getMessage());
                }

                executeAndWait(f_script, logger);
            });

            rank++;
        }
    }
}
