package lib.Tgene;

import lib.ExecutableStep;
import tfprio.TFPRIO;
import util.Configs.Config;

import java.io.*;
import java.util.*;

import static util.FileManagement.extend;
import static util.FileManagement.makeSureFileExists;

public class FilterTargetGenes extends ExecutableStep
{
    private final Config<File> d_input_tgene = TFPRIO.configs.tgene.fileStructure.d_groups;
    private final Config<File> d_input_tepic = TFPRIO.configs.tepic.fileStructure.d_postprocessing_output;

    private final Config<File> d_output = TFPRIO.configs.tgene.fileStructure.d_filteredTargetGenes;

    private final Config<String> s_tgene = TFPRIO.configs.tgene.fileStructure.s_groups_mergedGroups;
    private final Config<String> s_tepicRatiosDir =
            TFPRIO.configs.tepic.fileStructure.s_postprocessing_output_ratiosDir;

    @Override protected Set<Config<File>> getRequiredFileStructure()
    {
        return new HashSet<>(Arrays.asList(d_input_tgene, d_input_tepic));
    }

    @Override protected Set<Config<File>> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(d_output));
    }

    @Override protected Set<Config<?>> getRequiredConfigs()
    {
        return new HashSet<>(Arrays.asList(s_tgene, s_tepicRatiosDir));
    }

    @Override protected void execute()
    {
        for (Map.Entry<String, Set<String>> combinationEntry : TFPRIO.groupCombinationsToHms.entrySet())
        {
            String groupPairing = combinationEntry.getKey();
            for (String hm : combinationEntry.getValue())
            {
                executorService.submit(() ->
                {
                    File input_data_TGENE = extend(d_input_tgene.get(), groupPairing, hm, s_tgene.get());
                    File input_data_TEPIC = extend(d_input_tepic.get(), groupPairing, hm, s_tepicRatiosDir.get(),
                            groupPairing + ".txt");

                    HashSet<String> available_target_genes = new HashSet<>();

                    try (BufferedReader br_tgene = new BufferedReader(new FileReader(input_data_TGENE));)
                    {
                        String line_tgene;
                        br_tgene.readLine();
                        while ((line_tgene = br_tgene.readLine()) != null)
                        {
                            String[] split = line_tgene.split("\t");

                            String[] geneIDs = split[3].split(";");

                            for (String geneID : geneIDs)
                            {
                                if (TFPRIO.mapSymbolAndEnsg.hasEnsg(geneID))
                                {
                                    available_target_genes.add(geneID);
                                }
                            }
                        }
                    } catch (IOException e)
                    {
                        logger.error(e.getMessage());
                    }

                    File targetFile = extend(d_output.get(), groupPairing, hm, input_data_TEPIC.getName());
                    makeSureFileExists(targetFile, logger);

                    try (BufferedWriter bw_tepic = new BufferedWriter(new FileWriter(targetFile));
                         BufferedReader br_tepic = new BufferedReader(new FileReader(input_data_TEPIC)))
                    {
                        String line_tepic = br_tepic.readLine();

                        bw_tepic.write(line_tepic);
                        bw_tepic.newLine();

                        while ((line_tepic = br_tepic.readLine()) != null)
                        {
                            String geneID = line_tepic.substring(0, line_tepic.indexOf("\t"));
                            if (available_target_genes.contains(geneID))
                            {
                                bw_tepic.write(line_tepic);
                                bw_tepic.newLine();
                            }
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
