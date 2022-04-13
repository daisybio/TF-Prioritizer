package lib.Tepic;

import lib.ExecutableStep;
import tfprio.TFPRIO;
import util.Configs.Config;
import util.FileFilters.Filters;

import java.io.*;
import java.util.*;

import static util.FileManagement.extend;

public class Randomize extends ExecutableStep
{
    private final Config<File> d_input = TFPRIO.configs.tepic.fileStructure.d_outputRaw;
    private final Config<File> d_output = TFPRIO.configs.tepic.fileStructure.d_outputRaw_shuffle;
    private final Config<String> s_trapSequences = TFPRIO.configs.tepic.fileStructure.s_outputRaw_trapSequences;

    @Override protected Set<Config<File>> getRequiredFileStructure()
    {
        return new HashSet<>(List.of(d_input));
    }

    @Override protected Set<Config<File>> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(d_output));
    }

    @Override protected Set<Config<?>> getRequiredConfigs()
    {
        return new HashSet<>(List.of(s_trapSequences));
    }

    @Override protected void execute()
    {
        Random random = new Random();

        for (File d_group : Objects.requireNonNull(d_input.get().listFiles(Filters.directoryFilter)))
        {
            for (File d_hm : Objects.requireNonNull(d_group.listFiles(Filters.directoryFilter)))
            {
                for (File d_sample : Objects.requireNonNull(d_hm.listFiles(Filters.fileFilter)))
                {
                    executorService.execute(() ->
                    {
                        File d_outputSample =
                                extend(d_output.get(), d_group.getName(), d_hm.getName(), d_sample.getName());
                        File f_sample = extend(d_sample, s_trapSequences.get());

                        ArrayList<File> affinityFiles = new ArrayList<>();

                        for (File f_affinity : Objects.requireNonNull(d_sample.listFiles(Filters.fileFilter)))
                        {
                            if (f_affinity.getName().matches(".*_Affinity_Gene_View_Filtered.*"))
                            {
                                affinityFiles.add(f_affinity);
                            }
                        }

                        File file_to_shuffle = null;
                        if (affinityFiles.size() > 1)
                        {
                            for (File f : affinityFiles)
                            {
                                if (f.getName().matches(".*TPM.txt"))
                                {
                                    file_to_shuffle = f;
                                    break;
                                }
                            }
                        } else
                        {
                            file_to_shuffle = affinityFiles.get(0);
                        }

                        //now shuffle columns of this file and write line by line to new file
                        assert file_to_shuffle != null;
                        File f_out_shuffled = extend(d_outputSample, file_to_shuffle.getName());

                        try (BufferedWriter bw = new BufferedWriter(new FileWriter(f_out_shuffled));
                             BufferedReader br = new BufferedReader(new FileReader(file_to_shuffle)))
                        {
                            String line = br.readLine();
                            String[] split_header = line.split("\t");

                            bw.write(line);
                            bw.newLine();

                            ArrayList<Integer> which_column_writes_which = new ArrayList<>();
                            which_column_writes_which.add(0);

                            ArrayList<Integer> not_shuffled = new ArrayList<>();
                            for (int i = 1; i < split_header.length; i++)
                            {
                                not_shuffled.add(i);
                            }

                            Collections.shuffle(not_shuffled);
                            which_column_writes_which.addAll(not_shuffled);

                            while ((line = br.readLine()) != null)
                            {
                                String[] split = line.split("\t");

                                StringBuilder sb = new StringBuilder();

                                int c = 0;
                                for (Integer i : which_column_writes_which)
                                {
                                    sb.append(split[i]);
                                    if (c < which_column_writes_which.size() - 1)
                                    {
                                        sb.append("\t");
                                    }
                                    c++;
                                }

                                bw.write(sb.toString());
                                bw.newLine();
                            }
                        } catch (IOException e)
                        {
                            logger.error(e.getMessage());
                        }

                        //shuffle sequences
                        File f_out_shuffled_sequences = extend(d_outputSample, f_sample.getName());

                        HashSet<String> tfs = new HashSet<>();
                        try (BufferedReader br_sequences = new BufferedReader(new FileReader(f_sample)))
                        {
                            String line_shuf_seq;
                            while ((line_shuf_seq = br_sequences.readLine()) != null)
                            {
                                if (line_shuf_seq.startsWith("#"))
                                {
                                    continue;
                                }

                                String[] split = line_shuf_seq.split("\t");
                                if (tfs.contains(split[0]))
                                {
                                    break;
                                } else
                                {
                                    tfs.add(split[0]);
                                }
                            }
                        } catch (IOException e)
                        {
                            logger.error(e.getMessage());
                        }

                        ArrayList<String> tfs_draw_random = new ArrayList<>(tfs);

                        try (BufferedReader br_sequences_write = new BufferedReader(new FileReader(f_sample));
                             BufferedWriter bw_shuffled_seq = new BufferedWriter(
                                     new FileWriter(f_out_shuffled_sequences)))
                        {
                            String line_shuf_seq_write;
                            while ((line_shuf_seq_write = br_sequences_write.readLine()) != null)
                            {
                                if (line_shuf_seq_write.startsWith("#") || line_shuf_seq_write.startsWith("TF\t"))
                                {
                                    bw_shuffled_seq.write(line_shuf_seq_write);
                                    bw_shuffled_seq.newLine();
                                    continue;
                                }

                                String[] split = line_shuf_seq_write.split("\t");
                                StringBuilder sb_suf = new StringBuilder();

                                //draw random TF
                                int random_index = random.nextInt(tfs_draw_random.size() - 1);
                                String tf_drawn = tfs_draw_random.get(random_index);

                                sb_suf.append(tf_drawn);
                                for (int i = 1; i < split.length; i++)
                                {
                                    sb_suf.append("\n");
                                    sb_suf.append(split[i]);
                                }
                                bw_shuffled_seq.write(sb_suf.toString());
                                bw_shuffled_seq.newLine();
                            }
                        } catch (IOException e)
                        {
                            logger.error(e.getMessage());
                        }
                    });
                }
            }
        }
    }
}
