package lib.AngularReport;

import org.json.JSONObject;
import tfprio.TFPRIO;
import util.Logger;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

import static util.FileManagement.extend;

public class TranscriptionFactorGroup
{
    final String name;
    final Logger logger;

    final List<String> geneIDs = new ArrayList<>();

    final List<TranscriptionFactor> transcriptionFactors = new ArrayList<>();
    Set<TargetGene> targetGenes = new HashSet<>();

    public TranscriptionFactorGroup(String name, Logger logger)
    {
        this.name = name;
        this.logger = logger;

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

    public JSONObject toJSONObject()
    {
        return new JSONObject()
        {{
            accumulate("transcriptionFactors", new JSONObject()
            {{
                for (TranscriptionFactor transcriptionFactor : transcriptionFactors)
                {
                    accumulate(transcriptionFactor.name, transcriptionFactor.toJSONObject());
                }
            }});

            List<JSONObject> targetGeneJsonObjects = new ArrayList<>();
            targetGenes.forEach(targetGene -> targetGeneJsonObjects.add(targetGene.toJSONObject()));
            accumulate("targetGenes", targetGeneJsonObjects);
            accumulate("geneIDs", geneIDs);
        }};
    }
}
