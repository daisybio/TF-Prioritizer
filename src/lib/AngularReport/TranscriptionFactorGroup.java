package lib.AngularReport;

import org.json.JSONObject;
import tfprio.TFPRIO;
import util.FileFilters.Filters;
import util.Logger;

import javax.print.attribute.standard.JobStateReasons;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ExecutorService;

import static util.FileManagement.extend;
import static util.FileManagement.getFileIfInDirectory;

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
    JSONObject validation_logos_biophysicalModel;
    JSONObject validation_logos_tfSequence;

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

            Generate.linkFiles(validation_heatmap, extend(d_validation, "heatmap"), executorService, logger);
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

            Generate.linkFiles(validation_igv, extend(d_validation, "igv"), executorService, logger);
        } // IGV

        {
            {
                File d_source = TFPRIO.configs.distributionAnalysis.fileStructure.d_logos_biophysicalModel.get();

                File d_tfGroup = getFileIfInDirectory(d_source, "[0-9]+_" + name, false);

                Map<String, File> fileType_file = new HashMap<>();

                if (d_tfGroup != null)
                {
                    File f_plot = extend(d_tfGroup, "biophysical_model.png");
                    File f_data = extend(d_tfGroup, "energy_matrix.csv");

                    if (f_plot.exists())
                    {
                        fileType_file.put("plot", f_plot);
                    }
                    if (f_data.exists())
                    {
                        fileType_file.put("data", f_data);
                    }
                }

                validation_logos_biophysicalModel = new JSONObject(fileType_file);

                Generate.linkFiles(validation_logos_biophysicalModel, extend(d_validation, "logos", "biophysical"),
                        executorService, logger);
            } // Biophysical model

            {
                File d_source = TFPRIO.configs.distributionAnalysis.fileStructure.d_logos_tfSequence.get();

                File d_tfGroup = getFileIfInDirectory(d_source, "[0-9]+_" + name, false);

                Map<String, Map<String, File>> subtype_filetype_file = new HashMap<>();

                if (d_tfGroup != null)
                {
                    Set<String> subtypes = new HashSet<>();

                    for (File f_found : Objects.requireNonNull(d_tfGroup.listFiles(Filters.fileFilter)))
                    {
                        String subtype = f_found.getName().substring(0, f_found.getName().lastIndexOf("."));
                        subtypes.add(subtype);
                    }

                    for (String subtype : subtypes)
                    {
                        subtype_filetype_file.put(subtype, new HashMap<>());

                        File f_data = extend(d_tfGroup, subtype + ".json");
                        File f_plot = extend(d_tfGroup, subtype + ".svg");

                        if (f_plot.exists())
                        {
                            subtype_filetype_file.get(subtype).put("plot", f_plot);
                        }
                        if (f_data.exists())
                        {
                            subtype_filetype_file.get(subtype).put("data", f_data);
                        }
                    }
                }

                validation_logos_tfSequence = new JSONObject(subtype_filetype_file);

                Generate.linkFiles(validation_logos_tfSequence, extend(d_validation, "logos", "tfSequence"),
                        executorService, logger);
            } // TF Sequence (JASPAR)
        } // Logos
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
                put("logos", new JSONObject()
                {{
                    put("biophysical", validation_logos_biophysicalModel);
                    put("tfSequence", validation_logos_tfSequence);
                }});
            }});
            put("targetGenes", targetGeneJsonObjects);
        }};
    }
}
