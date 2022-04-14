package lib.Tgene;

import lib.ExecutableStep;
import lib.GeneRegion;
import tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.FileFilters.Filters;

import java.io.*;
import java.util.*;

import static lib.GeneRegion.removeGeneRegionDuplicates;
import static util.FileManagement.*;

public class Postprocessing extends ExecutableStep
{
    private final AbstractConfig<File> d_input = TFPRIO.configs.tgene.fileStructure.d_output;
    private final GeneratedFileStructure d_output = TFPRIO.configs.tgene.fileStructure.d_merged;
    private final GeneratedFileStructure d_outputGroups = TFPRIO.configs.tgene.fileStructure.d_groups;
    private final GeneratedFileStructure f_transcripts_gtf = TFPRIO.configs.tgene.fileStructure.f_transcripts_gtf;

    private final AbstractConfig<File> d_preprocessingMeanCounts =
            TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_meanCounts;

    private final AbstractConfig<String> s_output_links = TFPRIO.configs.tgene.fileStructure.s_output_links;
    private final AbstractConfig<String> s_preprocessing_meanCounts =
            TFPRIO.configs.deSeq2.fileStructure.s_preprocessing_MeanCounts;
    private final AbstractConfig<String> s_groups_mergedGroups =
            TFPRIO.configs.tgene.fileStructure.s_groups_mergedGroups;
    private final AbstractConfig<Integer> tpmCutoff = TFPRIO.configs.tepic.tpmCutoff;
    private final AbstractConfig<Boolean> mutuallyExclusive = TFPRIO.configs.mixOptions.mutuallyExclusive;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(Arrays.asList(d_input, d_preprocessingMeanCounts, f_transcripts_gtf));
    }

    @Override public Set<GeneratedFileStructure> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(d_output))
        {{
            if (mutuallyExclusive.get())
            {
                add(d_outputGroups);
            } else
            {
                d_outputGroups.deleteAndSetNoGenerationReason(mutuallyExclusive.getName() + " is set to false");
            }
        }};
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>(
                Arrays.asList(s_output_links, s_preprocessing_meanCounts, s_groups_mergedGroups, mutuallyExclusive));
    }

    @Override protected Set<AbstractConfig<?>> getOptionalConfigs()
    {
        return new HashSet<>(List.of(tpmCutoff));
    }

    @Override protected void execute()
    {
        Map<String, Integer> geneID_lengths = getGeneLengths();

        for (File d_group : Objects.requireNonNull(d_input.get().listFiles(Filters.directoryFilter)))
        {
            Map<String, Integer> geneID_counts = getGeneCounts(d_group.getName());
            double all_rpk = getRpkSum(geneID_counts, geneID_lengths);
            double scaling_factor = all_rpk / 1E6;

            for (File d_hm : Objects.requireNonNull(d_group.listFiles(Filters.directoryFilter)))
            {
                executorService.submit(() ->
                {
                    File d_output_groupHm = extend(d_output.get(), d_group.getName(), d_hm.getName());

                    List<GeneRegion> allRegions = new ArrayList<>();
                    Map<String, List<GeneRegion>> sampleSpecificRegions = new HashMap<>();

                    for (File d_sample : Objects.requireNonNull(d_hm.listFiles(Filters.directoryFilter)))
                    {
                        //this one is needed for mutually exclusive option
                        String sampleName = d_sample.getName();
                        List<GeneRegion> currentSampleRegions = new ArrayList<>();

                        File f_links = extend(d_sample, s_output_links.get());

                        try (BufferedReader br = new BufferedReader(new FileReader(f_links)))
                        {
                            String inputLine;
                            br.readLine();
                            int lineCount = 0;

                            while ((inputLine = br.readLine()) != null)
                            {
                                if (inputLine.startsWith("#") || inputLine.isBlank())
                                {
                                    continue;
                                }

                                String[] inputSplit = inputLine.split("\t");

                                String ensg = inputSplit[0].substring(0, inputSplit[0].indexOf("."));
                                String chromosome = inputSplit[6].substring(0, inputSplit[6].indexOf(":"));
                                int start = Integer.parseInt(inputSplit[6].substring(inputSplit[6].indexOf(":") + 1,
                                        inputSplit[6].indexOf("-")));
                                int end = Integer.parseInt(inputSplit[6].substring(inputSplit[6].indexOf("-") + 1));

                                GeneRegion region = new GeneRegion(chromosome, start, end, ensg);

                                currentSampleRegions.add(region);
                                lineCount++;
                            }
                        } catch (IOException e)
                        {
                            e.printStackTrace();
                            System.exit(1);
                        }
                        List<GeneRegion> currentSampleRegionsNoDuplicates =
                                removeGeneRegionDuplicates(currentSampleRegions, true);
                        allRegions.addAll(currentSampleRegionsNoDuplicates);
                        sampleSpecificRegions.put(sampleName, currentSampleRegionsNoDuplicates);
                    }
                    List<GeneRegion> allRegionsNoDuplicates = removeGeneRegionDuplicates(allRegions, true);
                    List<GeneRegion> allRegionsFiltered = allRegionsNoDuplicates;

                    if (tpmCutoff.isSet())
                    {
                        logger.info("TPM filter is set. Filtering TGENE results: " + d_hm.getName() + " - " +
                                d_group.getName() + ".");

                        allRegionsFiltered =
                                filterRegionsByRpkmThreshold(allRegionsNoDuplicates, scaling_factor, geneID_counts,
                                        geneID_lengths);

                        for (String sampleName : sampleSpecificRegions.keySet())
                        {
                            List<GeneRegion> currentSampleRegion = sampleSpecificRegions.get(sampleName);
                            List<GeneRegion> currentSampleRegionsFiltered =
                                    filterRegionsByRpkmThreshold(currentSampleRegion, scaling_factor, geneID_counts,
                                            geneID_lengths);
                            sampleSpecificRegions.put(sampleName, currentSampleRegionsFiltered);
                        }
                    }

                    if (!mutuallyExclusive.get())
                    {
                        String targetFileName = d_hm.getName() + "_" + d_group.getName() + ".tsv";
                        File targetFile = extend(d_output_groupHm, targetFileName);
                        writeGeneRegionsToFile(allRegionsFiltered, targetFile);
                    } else
                    {
                        for (String sampleName : sampleSpecificRegions.keySet())
                        {
                            List<GeneRegion> regions = sampleSpecificRegions.get(sampleName);

                            String targetFileName = d_hm.getName() + "_" + sampleName.split("_")[0] + ".tsv";
                            File targetFile = extend(d_output_groupHm, targetFileName);

                            writeGeneRegionsToFile(regions, targetFile);
                        }

                        //create tgen groups
                        File targetFile = extend(d_outputGroups.get(), d_hm.getName(), d_group.getName(),
                                s_groups_mergedGroups.get());
                        writeGeneRegionsToFile(allRegionsFiltered, targetFile);
                    }
                });
            }
        }
    }

    private List<GeneRegion> filterRegionsByRpkmThreshold(List<GeneRegion> regions, double threshold,
                                                          Map<String, Integer> counts, Map<String, Integer> lengths)
    {
        List<GeneRegion> currentSampleRegionsFiltered = new ArrayList<>();

        for (GeneRegion region : regions)
        {
            //check if we want that TF

            String geneID = region.getId();

            if (!TFPRIO.mapSymbolAndEnsg.hasEnsg(geneID) || !counts.containsKey(geneID) || !lengths.containsKey(geneID))
            {
                continue;
            }

            int count = counts.get(geneID);
            int length = lengths.get(geneID);
            if (length == 0)
            {
                continue;
            }
            double lengthPerK = length / 1000.0;
            double rpk = count / lengthPerK;

            if (rpk >= threshold)
            {
                currentSampleRegionsFiltered.add(region);
            }
        }
        return currentSampleRegionsFiltered;
    }

    private Map<String, Integer> getGeneCounts(String groups)
    {
        Map<String, Integer> counts = new HashMap<>();
        if (!mutuallyExclusive.get())
        {
            File f_geneCount = extend(d_preprocessingMeanCounts.get(), groups + s_preprocessing_meanCounts.get());
            try (BufferedReader reader = new BufferedReader(new FileReader(f_geneCount)))
            {
                String inputLine;
                reader.readLine();
                while ((inputLine = reader.readLine()) != null)
                {
                    String[] split = inputLine.split("\t");
                    counts.put(split[1], Integer.parseInt(split[2]));
                }
            } catch (IOException e)
            {
                e.printStackTrace();
                System.exit(1);
            }
        } else
        {
            String[] names = groups.split("_");

            File f_geneCountGroup1 =
                    extend(d_preprocessingMeanCounts.get(), names[0] + s_preprocessing_meanCounts.get());
            File f_geneCountGroup2 =
                    extend(d_preprocessingMeanCounts.get(), names[1] + s_preprocessing_meanCounts.get());

            try (BufferedReader reader_geneCounts1 = new BufferedReader(new FileReader(f_geneCountGroup1));
                 BufferedReader reader_geneCounts2 = new BufferedReader(new FileReader(f_geneCountGroup2)))
            {
                String inputLine1, inputLine2;
                reader_geneCounts1.readLine();
                reader_geneCounts2.readLine();

                while ((inputLine1 = reader_geneCounts1.readLine()) != null &&
                        (inputLine2 = reader_geneCounts2.readLine()) != null)
                {
                    String[] split1 = inputLine1.split("\t");
                    String[] split2 = inputLine2.split("\t");
                    String id = split1[1];
                    assert split2[1].equals(id);
                    int count1 = Integer.parseInt(split1[2]);
                    int count2 = Integer.parseInt(split2[2]);
                    counts.put(id, (count1 + count2) / 2);
                }
            } catch (IOException e)
            {
                e.printStackTrace();
                System.exit(1);
            }
        }

        return counts;
    }

    private void writeGeneRegionsToFile(List<GeneRegion> regions, File targetFile)
    {
        try
        {
            makeSureFileExists(targetFile);
        } catch (IOException e)
        {
            e.printStackTrace();
            System.exit(1);
        }

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(targetFile)))
        {
            writer.write("CHR\tLEFT_BORDER\tRIGHT_BORDER\tENSG");
            writer.newLine();

            for (GeneRegion region : regions)
            {
                writer.write(region.toString());
                writer.newLine();

            }
        } catch (IOException e)
        {
            e.printStackTrace();
            System.exit(1);
        }
    }

    private double getRpkSum(Map<String, Integer> counts, Map<String, Integer> lengths)
    {
        double sum = 0.0;

        for (String geneID : counts.keySet())
        {
            if (!lengths.containsKey(geneID))
            {
                continue;
            }

            int count = counts.get(geneID);
            int length = lengths.get(geneID);
            if (length == 0)
            {
                continue;
            }
            double normalizedLength = length / 1000.0;
            double rpk = count / normalizedLength;
            sum += rpk;
        }

        return sum;
    }

    private Map<String, Integer> getGeneLengths()
    {
        HashMap<String, Integer> geneID_lengths = new HashMap<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(f_transcripts_gtf.get())))
        {
            String line_gene_symbol_lengths;
            while ((line_gene_symbol_lengths = reader.readLine()) != null)
            {
                if (line_gene_symbol_lengths.startsWith("#"))
                {
                    continue;
                }
                String[] split = line_gene_symbol_lengths.split("\t");
                int start = Integer.parseInt(split[3]);
                int end = Integer.parseInt(split[4]);
                int length = end - start + 1;
                String geneID = split[8].substring(split[8].indexOf("\"") + 1, split[8].indexOf("."));
                geneID_lengths.put(geneID, length);
            }
        } catch (IOException e)
        {
            e.printStackTrace();
            System.exit(1);
        }

        return geneID_lengths;
    }
}
