package lib.MixOptions;

import tfprio.TFPRIO;
import lib.ExecutableStep;
import util.FileFilters.Filters;

import java.io.*;
import java.util.*;
import java.util.concurrent.Executors;

import static util.FileManagement.*;

public class MixOptions extends ExecutableStep
{
    /**
     * preprocess mix histones, search for same peaks and use either the union or the intersection of all
     */
    public void execute()
    {
        logger.info("Used data: " + TFPRIO.configs.mixOptions.fileStructure.d_preprocessingCheckChr.get());
        preprocess();
        mainStep();
        if (TFPRIO.configs.mixOptions.level.get().equals("HM_LEVEL"))
        {
            hmLevel();
        }
    }

    private void splitFileByChromosome(File sourceFile, File targetDir)
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
            e.printStackTrace();
            System.exit(1);
        }

    }

    private void preprocess()
    {
        logger.info("Preprocess input data for sample mix - split chromosomes.");

        for (File d_group : Objects.requireNonNull(TFPRIO.configs.mixOptions.fileStructure.d_preprocessingCheckChr.get()
                .listFiles(Filters.directoryFilter)))
        {
            String group = d_group.getName();
            File d_outputGroup = extend(TFPRIO.configs.mixOptions.fileStructure.d_sampleMixPreprocessing.get(), group);

            for (File d_hm : d_group.listFiles(Filters.directoryFilter))
            {
                File d_outputGroupHm = extend(d_outputGroup, d_hm.getName());

                for (File f_sample : d_hm.listFiles(Filters.fileFilter))
                {
                    executorService.execute(() ->
                    {
                        splitFileByChromosome(f_sample, d_outputGroupHm);
                    });
                }
            }
        }

        shutdown();
        executorService = Executors.newFixedThreadPool(TFPRIO.configs.general.threadLimit.get());
    }

    private void mainStep()
    {
        runMixOption(TFPRIO.configs.mixOptions.fileStructure.d_sampleMixPreprocessing.get(),
                TFPRIO.configs.mixOptions.fileStructure.d_sampleMix.get(), "SAMPLE_LEVEL");
    }

    private void hmLevelPreprocess()
    {
        logger.info("Preprocess sample unions for HM mix");
        logger.info("Identify possible groups with same histone modifications");

        HashMap<String, ArrayList<String>> groups_histoneModifications = new HashMap<>();
        HashMap<String, ArrayList<String>> deleted_groups = new HashMap<>();
        HashSet<String> available_histoneModifications = new HashSet<>();

        //identify possible timepoints
        for (File d_group : Objects.requireNonNull(
                TFPRIO.configs.mixOptions.fileStructure.d_sampleMix.get().listFiles(Filters.directoryFilter)))
        {
            String group = d_group.getName();
            ArrayList<String> hmList = new ArrayList<>();

            for (File d_hm : Objects.requireNonNull(d_group.listFiles(Filters.directoryFilter)))
            {
                hmList.add(d_hm.getName());
                available_histoneModifications.add(d_hm.getName());
            }
            groups_histoneModifications.put(group, hmList);
        }


        for (String group : groups_histoneModifications.keySet())
        {
            if (groups_histoneModifications.get(group).size() < available_histoneModifications.size())
            {
                deleted_groups.put(group, groups_histoneModifications.get(group));
                continue;
            }

            if (!groups_histoneModifications.get(group).containsAll(available_histoneModifications))
            {
                deleted_groups.put(group, groups_histoneModifications.get(group));
            }
        }

        deleted_groups.keySet().forEach(groups_histoneModifications::remove);

        StringBuilder sb_foundMixingGroups = new StringBuilder();

        sb_foundMixingGroups.append("Can perform complete mix for HMs (");

        available_histoneModifications.forEach(
                histoneModification -> sb_foundMixingGroups.append(histoneModification).append(" "));

        sb_foundMixingGroups.append(") in timepoints (");

        groups_histoneModifications.keySet().forEach(group -> sb_foundMixingGroups.append(group).append(" "));

        sb_foundMixingGroups.append("). Can perform part-mix or no-mix for timepoints (");

        deleted_groups.keySet().forEach(deletedGroup -> sb_foundMixingGroups.append(deletedGroup).append(" "));

        sb_foundMixingGroups.append(").");
        logger.info(sb_foundMixingGroups.toString());


        for (File d_group : Objects.requireNonNull(
                TFPRIO.configs.mixOptions.fileStructure.d_sampleMix.get().listFiles(Filters.directoryFilter)))
        {
            String group = d_group.getName();
            File d_output_hmPreprocessing =
                    extend(TFPRIO.configs.mixOptions.fileStructure.d_preprocessingHmMix.get(), group);

            for (File d_hm : Objects.requireNonNull(d_group.listFiles(Filters.directoryFilter)))
            {
                File d_output_hmPreprocessingMix = extend(d_output_hmPreprocessing, "MIX");

                for (File f_sample : Objects.requireNonNull(d_hm.listFiles(Filters.fileFilter)))
                {
                    splitFileByChromosome(f_sample, d_output_hmPreprocessingMix);
                }
            }
        }
    }

    private void runMixOption(File d_source, File d_target, String mixLevel)
    {
        logger.info("Create " + TFPRIO.configs.mixOptions.option.get() + " of " + mixLevel);

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
                                System.exit(1);
                            }
                        }


                        //unions sorting
                        Collections.sort(all_intervals);

                        Stack<MIX_Interval> stack_union = mergeIntervals(all_intervals);

                        ArrayList<MIX_Interval> chromosomeUnions = new ArrayList<>();

                        int min_occurence;
                        if (!TFPRIO.configs.mixOptions.occurrenceIntersection.isSet())
                        {
                            min_occurence = sample_number;
                        } else
                        {
                            min_occurence =
                                    Math.min(sample_number, TFPRIO.configs.mixOptions.occurrenceIntersection.get());
                        }

                        while (!stack_union.isEmpty())
                        {
                            MIX_Interval t = stack_union.pop();
                            t.calculate_mean(mixLevel);

                            if (TFPRIO.configs.mixOptions.option.get().equals("INTERSECTION"))
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
                            System.exit(1);
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
                            System.exit(1);
                        }
                    }
                });
            }
        }

        try
        {
            TFPRIO.configs.general.latestInputDirectory.setValue(d_target);
        } catch (IllegalAccessException e)
        {
            logger.error(e.getMessage());
        }
    }

    private void hmLevel()
    {
        hmLevelPreprocess();

        runMixOption(TFPRIO.configs.mixOptions.fileStructure.d_preprocessingHmMix.get(),
                TFPRIO.configs.mixOptions.fileStructure.d_hmMix.get(), "HM_LEVEL");
    }


    /**
     * mixes the samples of one folder into one file, based on mix_option (UNION or INTERSECTION)
     */
    private static Stack<MIX_Interval> mergeIntervals(ArrayList<MIX_Interval> interval)
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
}