package lib.FootprintIntervals;

import tfprio.TFPRIO;
import lib.ExecutableStep;
import tfprio.Workflow;
import util.Comparators.ChromosomeComparator;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.FileFilters.Filters;

import java.io.*;
import java.util.*;

import static util.FileManagement.extend;
import static util.FileManagement.makeSureFileExists;

public class CreateFootprintsBetweenPeaks extends ExecutableStep
{
    private AbstractConfig<File> d_input;
    private final GeneratedFileStructure d_output = TFPRIO.configs.mixOptions.fileStructure.d_footprintsBetweenPeaks;
    private final AbstractConfig<String> option = TFPRIO.configs.tepic.tfBindingSiteSearch;


    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(List.of(d_input));
    }

    @Override public Set<GeneratedFileStructure> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(d_output));
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>(List.of(option));
    }

    @Override protected void updateInputDirectory()
    {
        d_input = Workflow.getLatestInputDirectory();
    }

    public void execute()
    {
        logger.info("OPTION tfBindingSiteSearch=\"" + option.get() + "\" was set.");
        logger.info("Creating footprints");


        logger.info("Used data: " + d_input.get().getAbsolutePath());


        for (File d_group : Objects.requireNonNull(d_input.get().listFiles(Filters.directoryFilter)))
        {
            String group = d_group.getName();
            File d_outputGroup = extend(d_output.get(), group);

            for (File d_hm : Objects.requireNonNull(d_group.listFiles(Filters.directoryFilter)))
            {
                String hm = d_hm.getName();
                File d_outputGroupHm = new File(d_outputGroup.getAbsolutePath() + File.separator + hm);

                for (File f_sample : Objects.requireNonNull(d_hm.listFiles(Filters.fileFilter)))
                {
                    executorService.execute(() ->
                    {
                        File f_outputPeak = extend(d_outputGroupHm, f_sample.getName());

                        //read in all regions in array list per chromosome
                        HashMap<String, ArrayList<Footprint_Interval>> chromosomeRegions = new HashMap<>();

                        try
                        {
                            try (BufferedReader reader_inputChromosomeRegions = new BufferedReader(
                                    new FileReader(f_sample)))
                            {
                                String inputLine;
                                while ((inputLine = reader_inputChromosomeRegions.readLine()) != null)
                                {
                                    String chromosome = inputLine.substring(0, inputLine.indexOf("\t"));

                                    if (!chromosomeRegions.containsKey(chromosome))
                                    {
                                        chromosomeRegions.put(chromosome, new ArrayList<>());
                                    }
                                }
                            }

                            try (BufferedReader reader_inputRegions = new BufferedReader(new FileReader(f_sample)))
                            {
                                String inputLine;
                                while ((inputLine = reader_inputRegions.readLine()) != null)
                                {
                                    String[] split = inputLine.split("\t");

                                    String chromosome = split[0];
                                    int start = Integer.parseInt(split[1]);
                                    int end = Integer.parseInt(split[2]);
                                    String name = split[3];
                                    int score = Integer.parseInt(split[4]);
                                    String strand = split[5];
                                    double signalValue = Double.parseDouble(split[6]);
                                    double pValue = Double.parseDouble(split[7]);
                                    double qValue = Double.parseDouble(split[8]);
                                    Footprint_Interval fp =
                                            new Footprint_Interval(chromosome, start, end, name, score, strand,
                                                    signalValue, pValue, qValue, inputLine);

                                    chromosomeRegions.get(chromosome).add(fp);
                                }
                            }
                        } catch (IOException e)
                        {
                            logger.error("Error during read process: " + e.getMessage());
                            System.exit(1);
                        }

                        //sort region array lists per chromosome
                        for (String chromosome : chromosomeRegions.keySet())
                        {
                            Collections.sort(chromosomeRegions.get(chromosome));
                        }

                        //determine chromosome order
                        ArrayList<String> chromosomesOrdered = new ArrayList<>(chromosomeRegions.keySet());
                        chromosomesOrdered.sort(new ChromosomeComparator());


                        //see if gaps between Footprint Intervals are < tepic_between_max_bps if so connect to one region for TEPIC

                        if (option.get().equals("BETWEEN"))
                        {
                            for (String chromosome : chromosomeRegions.keySet())
                            {
                                ArrayList<Footprint_Interval> intervals = chromosomeRegions.get(chromosome);

                                ArrayList<Footprint_Interval> intervals_combined = new ArrayList<>();

                                for (int i = 0; i < intervals.size(); i++)
                                {
                                    int startIndex = i;
                                    int currentIndex = i;

                                    int firstStartPosition = intervals.get(i).start;
                                    int lastEndPosition = intervals.get(i).end;

                                    StringBuilder name = new StringBuilder(intervals.get(i).name);
                                    StringBuilder strand = new StringBuilder(intervals.get(i).strand);

                                    int scoreSum = intervals.get(i).score;
                                    double signalValueSum = intervals.get(i).signalValue;
                                    double pValueSum = intervals.get(i).pValue;
                                    double qValueSum = intervals.get(i).qValue;

                                    boolean include_next_peak = true;

                                    while (include_next_peak)
                                    {
                                        if (currentIndex + 1 < intervals.size())
                                        {
                                            Footprint_Interval next = intervals.get(currentIndex + 1);

                                            if (next.start - lastEndPosition < TFPRIO.configs.tepic.betweenMaxBps.get())
                                            {
                                                lastEndPosition = next.end;
                                                name.append(";").append(next.name);
                                                strand.append(";").append(next.strand);

                                                scoreSum += next.score;
                                                signalValueSum += next.signalValue;
                                                pValueSum += next.pValue;
                                                qValueSum += next.qValue;

                                                currentIndex += 1;
                                            } else
                                            {
                                                include_next_peak = false;
                                            }
                                        } else
                                        {
                                            include_next_peak = false;
                                        }
                                    }

                                    if (startIndex == currentIndex)
                                    {
                                        intervals_combined.add(intervals.get(i));
                                    } else
                                    {
                                        int divisor = currentIndex - startIndex + 1;

                                        scoreSum /= divisor;
                                        signalValueSum /= divisor;
                                        pValueSum /= divisor;
                                        qValueSum /= divisor;

                                        Footprint_Interval combined =
                                                new Footprint_Interval(chromosome, firstStartPosition, lastEndPosition,
                                                        name.toString(), scoreSum, strand.toString(), signalValueSum,
                                                        pValueSum, qValueSum);
                                        intervals_combined.add(combined);

                                        i = currentIndex;
                                    }
                                }
                                chromosomeRegions.put(chromosome, intervals_combined);
                            }
                        }

                        if (option.get().equals("EXCL_BETWEEN"))
                        {
                            for (String key_chr : chromosomeRegions.keySet())
                            {
                                ArrayList<Footprint_Interval> intervals = chromosomeRegions.get(key_chr);

                                ArrayList<Footprint_Interval> intervals_combined = new ArrayList<>();

                                for (int i = 0; i < intervals.size(); i++)
                                {
                                    Footprint_Interval first = intervals.get(i);

                                    String name = first.name;
                                    int score = first.score;
                                    String strand = first.strand;
                                    double signalValue = first.signalValue;
                                    double pValue = first.pValue;
                                    double qValue = first.qValue;

                                    if (i + 1 < intervals.size() && (intervals.get(i + 1)).start - first.end <
                                            TFPRIO.configs.tepic.betweenMaxBps.get())
                                    {
                                        Footprint_Interval next = intervals.get(i + 1);

                                        name += ";" + next.name;
                                        strand += ";" + next.strand;

                                        score += next.score;
                                        signalValue += next.signalValue;
                                        pValue += next.pValue;
                                        qValue += next.qValue;

                                        i++;

                                        score /= 2;
                                        signalValue /= 2;
                                        pValue /= 2;
                                        qValue /= 2;

                                        Footprint_Interval combined =
                                                new Footprint_Interval(key_chr, first.end, next.start, name, score,
                                                        strand, signalValue, pValue, qValue);

                                        intervals_combined.add(combined);
                                    } else
                                    {
                                        intervals_combined.add(first);
                                    }
                                }
                                chromosomeRegions.put(key_chr, intervals_combined);
                            }
                        }


                        //write output file
                        try
                        {
                            makeSureFileExists(f_outputPeak);
                            BufferedWriter writer = new BufferedWriter(new FileWriter(f_outputPeak));

                            for (String chromosome : chromosomesOrdered)
                            {
                                ArrayList<Footprint_Interval> intervals = chromosomeRegions.get(chromosome);

                                for (Footprint_Interval interval : intervals)
                                {
                                    writer.write(interval.toString());
                                    writer.newLine();
                                }
                            }
                            writer.close();
                        } catch (IOException e)
                        {
                            logger.error("Error during write process: " + e.getMessage());
                            System.exit(1);
                        }
                    });
                }
            }
        }
    }
}
