package lib.Tgene;

import lib.ExecutableStep;
import lib.GeneRegion;
import tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;

import java.io.*;
import java.util.*;

import static lib.GeneRegion.removeGeneRegionDuplicates;
import static util.FileManagement.*;

public class CreateGroups extends ExecutableStep
{
    private final AbstractConfig<File> d_input = TFPRIO.configs.tgene.fileStructure.d_merged;

    private final AbstractConfig<File> d_output = TFPRIO.configs.tgene.fileStructure.d_groups;

    private final AbstractConfig<String> s_tgeneOutputGroups = TFPRIO.configs.tgene.fileStructure.s_groups_mergedGroups;
    private final AbstractConfig<Boolean> mutuallyExclusive = TFPRIO.configs.mixOptions.mutuallyExclusive;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(List.of(d_input));
    }

    @Override public Set<AbstractConfig<File>> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(d_output));
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>(Arrays.asList(s_tgeneOutputGroups, mutuallyExclusive));
    }

    @Override protected void execute()
    {
        //based on structure build necessary files
        for (Map.Entry<String, Set<String>> pairingEntry : TFPRIO.groupCombinationsToHms.entrySet())
        {
            String pairing = pairingEntry.getKey();

            String[] split_name = pairing.split("_");
            String group1 = split_name[0];
            String group2 = split_name[1];

            for (String hm : pairingEntry.getValue())
            {
                File targetFile = extend(d_output.get(), pairing, hm, s_tgeneOutputGroups.get());

                executorService.submit(() ->
                {
                    makeSureFileExists(targetFile, logger);

                    try (BufferedWriter writer = new BufferedWriter(new FileWriter(targetFile)))
                    {
                        writer.write("CHR\tLEFT_BORDER\tRIGHT_BORDER\tENSGS");
                        writer.newLine();

                        File input_group1;
                        File input_group2;

                        if (!mutuallyExclusive.get())
                        {
                            input_group1 = extend(d_input.get(), group1, hm, hm + "_" + group1 + ".tsv");
                            input_group2 = extend(d_input.get(), group2, hm, hm + "_" + group2 + ".tsv");
                        } else
                        {
                            input_group1 = extend(d_input.get(), pairing, hm, hm + "_" + group1 + ".tsv");
                            input_group2 = extend(d_input.get(), pairing, hm, hm + "_" + group2 + ".tsv");
                        }

                        List<GeneRegion> unmergedRegions = new ArrayList<>();
                        getGeneRegions(input_group1, unmergedRegions);
                        getGeneRegions(input_group2, unmergedRegions);

                        List<GeneRegion> noDuplicates = removeGeneRegionDuplicates(unmergedRegions, false);

                        for (GeneRegion region : noDuplicates)
                        {
                            writer.write(region.toString());
                            writer.newLine();
                        }
                    } catch (IOException e)
                    {
                        logger.error(e.getMessage());
                    }
                });
            }
        }
    }

    private void getGeneRegions(File inputFile, List<GeneRegion> unmergedRegions)
    {
        try (BufferedReader reader = new BufferedReader(new FileReader(inputFile)))
        {
            String inputLine;
            reader.readLine();

            while ((inputLine = reader.readLine()) != null)
            {
                String[] split = inputLine.split("\t");
                String chr = split[0];
                int left_border = Integer.parseInt(split[1]);
                int right_border = Integer.parseInt(split[2]);
                String geneID = split[3];

                unmergedRegions.add(new GeneRegion(chr, left_border, right_border, geneID));
            }
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }
    }
}
