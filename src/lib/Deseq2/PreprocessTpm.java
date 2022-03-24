package lib.Deseq2;

import lib.ExecutableStep;
import tfprio.TFPRIO;
import util.Configs.Config;
import util.FileFilters.Filters;

import java.io.*;
import java.util.*;

import static util.FileManagement.extend;

public class PreprocessTpm extends ExecutableStep
{
    private final Config<File> d_input = TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_tpm_tpmResults;
    private final Config<File> d_combined = TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_combined;
    private final Config<File> f_output_copy = TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_combinedOriginal;

    private final Config<Double> tpmFilter = TFPRIO.configs.deSeq2.tpmFilter;

    @Override protected Set<Config<File>> getRequiredFileStructure()
    {
        return new HashSet<>(Arrays.asList(d_combined, d_input));
    }

    @Override protected Set<Config<File>> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(f_output_copy));
    }

    @Override protected Set<Config<?>> getRequiredConfigs()
    {
        return new HashSet<>(List.of(tpmFilter));
    }

    @Override protected void execute()
    {
        {
            logger.info("Filter genes with TPM under cutoff.");

            //read in TPMs
            HashMap<String, HashMap<String, Double>> tp_gene_tpm_value = new HashMap<>();

            for (File f_group : Objects.requireNonNull(d_input.get().listFiles(Filters.fileFilter)))
            {
                String group = f_group.getName().substring(0, f_group.getName().lastIndexOf("."));

                if (group.split("_").length > 1)
                {
                    String[] split = group.split("_");
                    group = split[0].substring(0, 1).toUpperCase() + split[0].substring(1) +
                            split[1].substring(0, 1).toUpperCase() + split[1].substring(1);
                }

                HashMap<String, Double> gene_tpm = new HashMap<>();

                try (BufferedReader reader = new BufferedReader(new FileReader(f_group)))
                {
                    String line = reader.readLine();
                    List<String> split_header = List.of(line.split("\t"));
                    int tpmIndex = split_header.indexOf("tpm");
                    while ((line = reader.readLine()) != null)
                    {
                        String[] split = line.split("\t");
                        if (split[tpmIndex].equals("NA"))
                        {
                            continue;
                        }
                        gene_tpm.put(split[0], Double.parseDouble(split[tpmIndex]));
                    }
                } catch (IOException e)
                {
                    e.printStackTrace();
                }
                tp_gene_tpm_value.put(group, gene_tpm);
            }

            for (File f_input : Objects.requireNonNull(d_combined.get().listFiles(Filters.directoryFilter)))
            {
                String[] groups = f_input.getName().split("_");
                String name1 = groups[0];
                String name2 = groups[1];

                if (groups.length > 2)
                {
                    int found_non = -1;

                    for (int i = 0; i < groups.length; i++)
                    {
                        if (groups[i].equals("non") || groups[i].equals("Non"))
                        {
                            found_non = i;
                        }
                    }

                    name1 = groups[0].substring(0, 1).toUpperCase() + groups[0].substring(1);
                    name2 = groups[1].substring(0, 1).toUpperCase() + groups[1].substring(1);

                    if (found_non == 0)
                    {
                        name1 = groups[0].substring(0, 1).toUpperCase() + groups[0].substring(1) +
                                groups[1].substring(0, 1).toUpperCase() + groups[1].substring(1);
                        name2 = groups[2].substring(0, 1).toUpperCase() + groups[2].substring(1);
                    } else if (found_non == 1)
                    {
                        name1 = groups[0].substring(0, 1).toUpperCase() + groups[0].substring(1);
                        name2 = groups[1].substring(0, 1).toUpperCase() + groups[1].substring(1) +
                                groups[2].substring(0, 1).toUpperCase() + groups[2].substring(1);
                    }

                }


                File f_output_copy_group = extend(f_output_copy.get(), f_input.getName());

                for (File f_input_read : Objects.requireNonNull(f_input.listFiles(Filters.fileFilter)))
                {
                    StringBuilder sb_copy = new StringBuilder();

                    StringBuilder sb_tpm_filtered = new StringBuilder();

                    try (BufferedReader br = new BufferedReader(new FileReader(f_input_read)))
                    {
                        String line = br.readLine();

                        sb_copy.append(line);
                        sb_copy.append("\n");

                        sb_tpm_filtered.append(line);
                        sb_tpm_filtered.append("\n");

                        String[] header = line.split("\t");

                        for (int i = 1; i < header.length; i++)
                        {
                            if (header[i].matches(name1 + ".*"))
                            {
                                header[i] = name1;
                            }
                            if (header[i].matches(name2 + ".*"))
                            {
                                header[i] = name2;
                            }
                        }

                        while ((line = br.readLine()) != null)
                        {
                            sb_copy.append(line);
                            sb_copy.append("\n");

                            String[] split = line.split("\t");

                            sb_tpm_filtered.append(split[0]);

                            boolean one_changed = false;

                            StringBuilder changer = new StringBuilder();

                            for (int i = 1; i < split.length; i++)
                            {
                                String group_name = header[i];
                                HashMap<String, Double> lookup = new HashMap<>();

                                if (tp_gene_tpm_value.containsKey(group_name))
                                {
                                    lookup = tp_gene_tpm_value.get(group_name);
                                }


                                if (!lookup.containsKey(split[0]))
                                {
                                    changer.append("\t");
                                    changer.append("0");
                                    one_changed = true;
                                    continue;
                                }

                                double tpm_value = lookup.get(split[0]);
                                if (tpm_value < tpmFilter.get())
                                {
                                    changer.append("\t");
                                    changer.append("0");
                                    one_changed = true;
                                } else
                                {
                                    changer.append("\t");
                                    changer.append(split[i]);
                                }
                            }
                            changer.append("\n");

                            if (one_changed)
                            {
                                String sb_complete_changer = "\t0".repeat(split.length - 1) + "\n";
                                sb_tpm_filtered.append(sb_complete_changer);

                            } else
                            {
                                sb_tpm_filtered.append(changer);
                            }

                        }
                    } catch (IOException e)
                    {
                        e.printStackTrace();
                    }

                    try (BufferedWriter bw_copy = new BufferedWriter(
                            new FileWriter(extend(f_output_copy_group, f_input_read.getName())));
                         BufferedWriter bw_tpm_filtered = new BufferedWriter(new FileWriter(f_input_read)))
                    {
                        bw_copy.write(sb_copy.toString());
                        bw_tpm_filtered.write(sb_tpm_filtered.toString());
                    } catch (IOException e)
                    {
                        e.printStackTrace();
                    }
                }
            }
        }
    }
}
