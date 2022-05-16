package lib.MixOptions;

import util.Configs.ConfigTypes.AbstractConfig;
import util.FileFilters.Filters;
import util.Logger;

import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutorService;

import static util.FileManagement.*;

public class StaticMethods
{
    static void runMixOption(File d_source, File d_target, String mixLevel, Logger logger,
                             ExecutorService executorService, AbstractConfig<String> option,
                             AbstractConfig<Integer> occurrenceIntersection)
    {
        logger.info("Create " + option.get() + " of " + mixLevel);

        for (File d_group : Objects.requireNonNull(d_source.listFiles(Filters.directoryFilter)))
        {
            String group = d_group.getName();
            File d_outputGroup = extend(d_target, group);

            for (File d_hm : Objects.requireNonNull(d_group.listFiles(Filters.directoryFilter)))
            {
                executorService.execute(() ->
                {
                    String hm = d_hm.getName();
                    File f_output_union_samples_tp_hm = extend(d_outputGroup, hm);

                    HashMap<String, ArrayList<MIX_Interval>> allChromosomes = new HashMap<>();
                    ArrayList<Integer> chromosomeAlpha = new ArrayList<>();
                    ArrayList<String> chromosomeStr = new ArrayList<>();

                    for (File d_chromosome : Objects.requireNonNull(d_hm.listFiles(Filters.directoryFilter)))
                    {
                        String chromosomeName = d_chromosome.getName();

                        ArrayList<MIX_Interval> all_intervals = new ArrayList<>();

                        int sample_number = 0;
                        String fileExtension = "";


                        for (File f_sample : Objects.requireNonNull(d_chromosome.listFiles(Filters.fileFilter)))
                        {
                            fileExtension = f_sample.getName().substring(f_sample.getName().lastIndexOf("."));

                            //build array of a union of all peaks

                            try (BufferedReader reader = new BufferedReader(new FileReader(f_sample)))
                            {
                                String inputLine;
                                while ((inputLine = reader.readLine()) != null)
                                {
                                    String[] split = inputLine.split("\t");

                                    MIX_Interval_Object mio =
                                            new MIX_Interval_Object(split[0], Integer.parseInt(split[1]),
                                                    Integer.parseInt(split[2]), split[3], Integer.parseInt(split[4]),
                                                    split[5], Double.parseDouble(split[6]),
                                                    Double.parseDouble(split[7]), Double.parseDouble(split[8]));
                                    MIX_Interval mi =
                                            new MIX_Interval(Integer.parseInt(split[1]), Integer.parseInt(split[2]));

                                    mi.merged_intervals.add(mio);
                                    all_intervals.add(mi);
                                }
                                sample_number++;
                            } catch (IOException e)
                            {
                                logger.error(e.getMessage());
                            }
                        }


                        //unions sorting
                        Collections.sort(all_intervals);

                        Stack<MIX_Interval> stack_union = mergeIntervals(all_intervals);

                        ArrayList<MIX_Interval> chromosomeUnions = new ArrayList<>();

                        int min_occurence;
                        if (!occurrenceIntersection.isSet())
                        {
                            min_occurence = sample_number;
                        } else
                        {
                            min_occurence = Math.min(sample_number, occurrenceIntersection.get());
                        }

                        while (!stack_union.isEmpty())
                        {
                            MIX_Interval t = stack_union.pop();
                            t.calculate_mean(mixLevel);

                            if (option.get().equals("INTERSECTION"))
                            {
                                if (t.merged_intervals.size() >= min_occurence)
                                {
                                    chromosomeUnions.add(t);
                                }
                            } else
                            {
                                chromosomeUnions.add(t);
                            }
                        }

                        Collections.sort(chromosomeUnions);

                        allChromosomes.put(chromosomeName, chromosomeUnions);

                        try
                        {
                            chromosomeAlpha.add(Integer.parseInt(chromosomeName));
                        } catch (Exception e)
                        {
                            chromosomeStr.add(chromosomeName);
                        }

                        Collections.sort(chromosomeAlpha);
                        Collections.sort(chromosomeStr);


                        //print Sample_unions => can be used if not HM_Level is used
                        File targetFile = extend(f_output_union_samples_tp_hm, group + "_" + hm + fileExtension);

                        try
                        {
                            makeSureFileExists(targetFile);
                        } catch (IOException e)
                        {
                            logger.error(e.getMessage());
                        }

                        try (BufferedWriter writer = new BufferedWriter(new FileWriter(targetFile)))
                        {
                            int peak_counter = 1;

                            for (int chr : chromosomeAlpha)
                            {
                                ArrayList<MIX_Interval> x = allChromosomes.get("" + chr);

                                for (MIX_Interval mix_interval : x)
                                {
                                    writer.write(mix_interval.meanToString(peak_counter));
                                    writer.newLine();
                                    peak_counter++;
                                }
                            }

                            for (String chr : chromosomeStr)
                            {
                                ArrayList<MIX_Interval> x = allChromosomes.get(chr);
                                for (MIX_Interval mix_interval : x)
                                {
                                    writer.write(mix_interval.meanToString(peak_counter));
                                    writer.newLine();
                                    peak_counter++;
                                }
                            }
                        } catch (IOException e)
                        {
                            logger.error("Could not write to file: " + targetFile.getAbsolutePath());
                        }
                    }
                });
            }
        }
    }

    /**
     * mixes the samples of one folder into one file, based on mix_option (UNION or INTERSECTION)
     */
    static Stack<MIX_Interval> mergeIntervals(ArrayList<MIX_Interval> interval)
    {
        Stack<MIX_Interval> stack = new Stack<>();

        if (stack.empty())
        {
            stack.push(interval.get(0));
        }

        for (int i = 1; i < interval.size(); i++)
        {
            MIX_Interval top = stack.peek();

            if (top.end < interval.get(i).start)
            {
                stack.push(interval.get(i));
            } else if (top.end < interval.get(i).end)
            {
                top.end = interval.get(i).end;
                top.merged_intervals.addAll(interval.get(i).merged_intervals);

                stack.pop();
                stack.push(top);
            } else
            {
                top.merged_intervals.addAll(interval.get(i).merged_intervals);

                stack.pop();
                stack.push(top);
            }
        }
        return stack;
    }

    static void splitFileByChromosome(File sourceFile, File targetDir, Logger logger)
    {
        BufferedWriter writer = null;
        String currentChromosome = null;
        boolean first = true;

        try (BufferedReader reader = new BufferedReader(new FileReader(sourceFile)))
        {
            String inputLine;

            while ((inputLine = reader.readLine()) != null)
            {
                String chromosome = inputLine.substring(0, inputLine.indexOf("\t"));

                if (first)
                {
                    currentChromosome = chromosome;
                    File targetFile = extend(extend(targetDir, currentChromosome), sourceFile.getName());
                    makeSureFileExists(targetFile);
                    writer = new BufferedWriter(new FileWriter(targetFile));
                    first = false;
                }

                if (chromosome.equals(currentChromosome))
                {
                    writer.write(inputLine);
                    writer.newLine();
                }
                if (!chromosome.equals(currentChromosome))
                {
                    writer.close();
                    currentChromosome = chromosome;
                    File targetFile = extend(extend(targetDir, currentChromosome), sourceFile.getName());
                    makeSureFileExists(targetFile);
                    writer = new BufferedWriter(new FileWriter(targetFile));
                    writer.write(inputLine);
                    writer.newLine();
                }
            }
            assert writer != null;
            writer.close();
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

    }
}