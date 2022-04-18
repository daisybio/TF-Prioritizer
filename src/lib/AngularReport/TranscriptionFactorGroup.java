package lib.AngularReport;

import org.json.JSONObject;
import tfprio.TFPRIO;
import util.FileFilters.Filters;
import util.Logger;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ExecutorService;

import static util.FileManagement.extend;

public class TranscriptionFactorGroup
{
    final String name;
    final Logger logger;
    final ExecutorService executorService;

    final List<String> geneIDs = new ArrayList<>();

    final List<TranscriptionFactor> transcriptionFactors = new ArrayList<>();
    Set<TargetGene> targetGenes = new HashSet<>();
    JSONObject validation_heatmap;
    JSONObject validation_igv;

    final File d_tfData;

    public TranscriptionFactorGroup(String name, Logger logger, ExecutorService executorService)
    {
        this.name = name;
        this.logger = logger;
        this.executorService = executorService;

        d_tfData = extend(TFPRIO.configs.angularReport.fileStructure.d_data.get(), name);

        for (String transcriptionFactorName : name.split("\\.\\."))
        {
            transcriptionFactors.add(new TranscriptionFactor(transcriptionFactorName, logger));
        }

        try
        {
            geneIDs.addAll(TFPRIO.mapSymbolAndEnsg.symbolToEnsg(name));
        } catch (NoSuchFieldException e)
        {
            logger.warn("No geneIDs found for " + name);
        }
    }

    public void collectData()
    {
        transcriptionFactors.forEach(TranscriptionFactor::collectData);
        collectTargetGenes();
        collectValidationFiles();
    }

    private void collectTargetGenes()
    {
        for (Map.Entry<String, Set<String>> entry : TFPRIO.groupCombinationsToHms.entrySet())
        {
            String groupCombination = entry.getKey();

            for (String hm : entry.getValue())
            {
                File f_input = extend(TFPRIO.configs.distributionAnalysis.fileStructure.d_heatmaps.get(), name, hm,
                        groupCombination + ".csv");

                Integer c_geneID = null, c_geneSymbol = null;

                try (BufferedReader reader = new BufferedReader(new FileReader(f_input)))
                {
                    String[] headerSplit = reader.readLine().split(",");

                    for (int i = 0; i < headerSplit.length; i++)
                    {
                        if (headerSplit[i].equals("\"Geneid\""))
                        {
                            c_geneID = i;
                        } else if (headerSplit[i].equals("\"geneSymbol\""))
                        {
                            c_geneSymbol = i;
                        }
                    }

                    if (c_geneID == null || c_geneSymbol == null)
                    {
                        logger.error("Could not find geneID or geneSymbol column in file " + f_input.getAbsolutePath());
                    } else
                    {
                        String inputLine;

                        while ((inputLine = reader.readLine()) != null)
                        {
                            String[] split = inputLine.split(",");
                            targetGenes.add(new TargetGene(split[c_geneID].replace("\"", ""),
                                    split[c_geneSymbol].replace("\"", "")));
                        }
                    }
                } catch (IOException e)
                {
                    logger.error(e.getMessage());
                }
            }
        }
    }

    private void collectValidationFiles()
    {
        File d_validation = extend(d_tfData, "validation");

        {
            Map<String, Map<String, Map<String, File>>> hm_groupPairing_filetype_file = new HashMap<>();
            File d_input = extend(TFPRIO.configs.distributionAnalysis.fileStructure.d_heatmaps.get(), name);

            for (String hm : TFPRIO.existingHms)
            {
                hm_groupPairing_filetype_file.put(hm, new HashMap<>());

                for (String groupPairing : TFPRIO.groupCombinationsToHms.keySet())
                {
                    hm_groupPairing_filetype_file.get(hm).put(groupPairing, new HashMap<>());

                    File f_data = extend(d_input, hm, groupPairing + ".csv");
                    File f_plot = extend(d_input, hm, groupPairing + ".png");

                    if (f_data.exists() && f_plot.exists())
                    {
                        hm_groupPairing_filetype_file.get(hm).get(groupPairing).put("data", f_data);
                        hm_groupPairing_filetype_file.get(hm).get(groupPairing).put("plot", f_plot);
                    }
                }
            }

            validation_heatmap = new JSONObject(hm_groupPairing_filetype_file);

            try
            {
                Generate.linkFiles(validation_heatmap, extend(d_validation, "heatmap"), executorService, logger);
            } catch (IOException e)
            {
                logger.error(e.getMessage());
            }
        } // Heatmaps

        {
            Map<String, Map<String, Map<String, Map<String, File>>>> hm_groupPairing_targetGene_filetype_file =
                    new HashMap<>();
            File d_input = extend(TFPRIO.configs.igv.fileStructure.d_igvDcgTargetGenes.get(), name);

            for (String hm : TFPRIO.existingHms)
            {
                hm_groupPairing_targetGene_filetype_file.put(hm, new HashMap<>());

                for (String groupPairing : TFPRIO.groupCombinationsToHms.keySet())
                {
                    hm_groupPairing_targetGene_filetype_file.get(hm).put(groupPairing, new HashMap<>());

                    File d_groupPairing = extend(d_input, hm, groupPairing);

                    if (d_groupPairing.exists())
                    {
                        for (File plotFile : d_groupPairing.listFiles(Filters.getSuffixFilter(".png")))
                        {
                            String targetGene = plotFile.getName().replaceAll("\\d+_", "").replace(".png", "");

                            hm_groupPairing_targetGene_filetype_file.get(hm).get(groupPairing)
                                    .put(targetGene, new HashMap<>()
                                    {{
                                        put("plot", plotFile);
                                    }});
                        }
                    }
                }
            }

            validation_igv = new JSONObject(hm_groupPairing_targetGene_filetype_file);

            try
            {
                Generate.linkFiles(validation_igv, extend(d_validation, "igv"), executorService, logger);
            } catch (IOException e)
            {
                logger.error(e.getMessage());
            }
        } // IGV
    }

    public JSONObject toJSONObject()
    {
        return new JSONObject()
        {{
            put("name", name);
            put("geneIDs", geneIDs);

            List<JSONObject> transcriptionFactorObjects = new ArrayList<>();

            for (TranscriptionFactor transcriptionFactor : transcriptionFactors)
            {
                transcriptionFactorObjects.add(transcriptionFactor.toJSONObject());
            }

            put("transcriptionFactors", transcriptionFactorObjects);

            List<JSONObject> targetGeneJsonObjects = new ArrayList<>();
            targetGenes.forEach(targetGene -> targetGeneJsonObjects.add(targetGene.toJSONObject()));

            put("validation", new JSONObject()
            {{
                put("heatmap", validation_heatmap);
                put("igv", validation_igv);
            }});
            put("targetGenes", targetGeneJsonObjects);
        }};
    }
}
