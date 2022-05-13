package lib.Tgene;

import lib.ExecutableStep;
import tfprio.tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;

import java.io.*;
import java.util.*;

import static util.FileManagement.extend;
import static util.FileManagement.makeSureFileExists;

public class FilterTargetGenes extends ExecutableStep
{
    private final AbstractConfig<File> d_input_tgene = TFPRIO.configs.tgene.fileStructure.d_groups;
    private final AbstractConfig<File> d_input_tepic = TFPRIO.configs.tepic.fileStructure.d_postprocessing_output;

    private final GeneratedFileStructure d_output = TFPRIO.configs.tgene.fileStructure.d_filteredTargetGenes;

    private final AbstractConfig<String> s_tgene = TFPRIO.configs.tgene.fileStructure.s_groups_mergedGroups;
    private final AbstractConfig<String> s_tepicRatiosDir =
            TFPRIO.configs.tepic.fileStructure.s_postprocessing_output_ratiosDir;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(Arrays.asList(d_input_tgene, d_input_tepic));
    }

    @Override public Set<GeneratedFileStructure> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(d_output));
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
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
                    File f_input_tgene = extend(d_input_tgene.get(), groupPairing, hm, s_tgene.get());
                    File f_input_tepic = extend(d_input_tepic.get(), groupPairing, hm, s_tepicRatiosDir.get(),
                            groupPairing + ".txt");

                    HashSet<String> available_target_genes = new HashSet<>();

                    try (BufferedReader br_tgene = new BufferedReader(new FileReader(f_input_tgene));)
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

                    File targetFile = extend(d_output.get(), groupPairing, hm, f_input_tepic.getName());
                    makeSureFileExists(targetFile, logger);

                    try (BufferedWriter writer = new BufferedWriter(new FileWriter(targetFile));
                         BufferedReader reader = new BufferedReader(new FileReader(f_input_tepic)))
                    {
                        String line_tepic = reader.readLine();

                        writer.write(line_tepic);
                        writer.newLine();

                        while ((line_tepic = reader.readLine()) != null)
                        {
                            String geneID = line_tepic.substring(0, line_tepic.indexOf("\t"));
                            if (available_target_genes.contains(geneID))
                            {
                                writer.write(line_tepic);
                                writer.newLine();
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
