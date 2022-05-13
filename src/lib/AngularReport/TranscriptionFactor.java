package lib.AngularReport;

import org.json.JSONObject;
import tfprio.tfprio.TFPRIO;
import util.FileManagement;
import util.Logger;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

import static util.FileManagement.extend;
import static util.FileManagement.findValueInTable;

public class TranscriptionFactor
{
    final String name;
    private final Logger logger;

    final List<String> separateTranscriptionFactors;
    final String geneID;

    Map<String, Map<String, Double>> log2fc = new HashMap<>();
    Map<String, Double> tpm = new HashMap<>();
    Map<String, Integer> normalizedExpression = new HashMap<>();

    public TranscriptionFactor(String name, Logger logger)
    {
        this.name = name;
        this.logger = logger;

        separateTranscriptionFactors = List.of(name.split(".."));
        Set<String> geneIDs = new HashSet<>();

        try
        {
            geneIDs.addAll(TFPRIO.mapSymbolAndEnsg.symbolToEnsg(name));
        } catch (NoSuchFieldException e)
        {
            logger.warn(e.getMessage());
        }

        if (geneIDs.size() == 0)
        {
            logger.warn("Could not find geneID for " + name + ".");

            geneID = "Not found";
        } else
        {
            geneID = new ArrayList<>(geneIDs).get(0);

            if (geneIDs.size() > 1)
            {
                logger.warn("Found more than one geneID for " + name);
            }
        }
    }

    public void collectData()
    {
        collectLog2fc();
        collectTpm();
        collectNormalizedExpression();
    }

    private void collectNormalizedExpression()
    {
        for (String group : TFPRIO.groupsToHms.keySet())
        {
            File f_input = extend(TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_meanCounts.get(), group + ".tsv");

            Integer value = null;

            try
            {
                value = Integer.parseInt(FileManagement.findValueInTable(geneID, 1, 2, f_input, "\t", false));
            } catch (FileNotFoundException e)
            {
                logger.error(e.getMessage());
            } catch (NoSuchFieldException e)
            {
                logger.warn(e.getMessage());
            }

            normalizedExpression.put(group, value);
        }
    }

    private void collectTpm()
    {
        for (String group : TFPRIO.groupsToHms.keySet())
        {
            File f_input =
                    extend(TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_tpm_tpmResults.get(), group + ".tsv");

            Double value = null;
            try
            {
                value = Double.parseDouble(FileManagement.findValueInTable(geneID, 1, 4, f_input, "\t", true));
            } catch (FileNotFoundException e)
            {
                logger.error(e.getMessage());
            } catch (NoSuchFieldException e)
            {
                logger.warn(e.getMessage());
            }
            tpm.put(group, value);
        }
    }

    private void collectLog2fc()
    {
        for (String group : TFPRIO.groupsToHms.keySet())
        {
            log2fc.put(group, new HashMap<>());
        }

        for (String group1 : TFPRIO.groupsToHms.keySet())
        {
            for (String group2 : TFPRIO.groupsToHms.keySet())
            {
                if (group2.compareTo(group1) > 0)
                {
                    String groupPairing = group1 + "_" + group2;
                    File f_log2fc = extend(TFPRIO.configs.deSeq2.fileStructure.d_output.get(), groupPairing + ".tsv");

                    Double value = null;

                    try
                    {
                        value = Double.parseDouble(findValueInTable(geneID, 0, 1, f_log2fc, "\t", true));
                    } catch (FileNotFoundException e)
                    {
                        logger.error(e.getMessage());
                    } catch (NoSuchFieldException ignore)
                    {
                    }

                    log2fc.get(group1).put(group2, value);
                }
            }
        }
    }

    public JSONObject toJSONObject()
    {
        return new JSONObject()
        {{
            put("name", name);
            put("geneID", geneID);

            put("log2fc", new ArrayList<>()
            {{
                for (String group1 : log2fc.keySet())
                {
                    for (String group2 : log2fc.get(group1).keySet())
                    {
                        add(new JSONObject()
                        {{
                            put("groups", Arrays.asList(group1, group2));
                            put("value", log2fc.get(group1).get(group2));
                        }});
                    }
                }
            }});


            put("tpm", new ArrayList<>()
            {{
                for (String group : tpm.keySet())
                {
                    add(new JSONObject()
                    {{
                        put("groups", List.of(group));
                        put("value", tpm.get(group));
                    }});
                }
            }});
            put("normalizedExpression", new ArrayList<>()
            {{
                for (String group : normalizedExpression.keySet())
                {
                    add(new JSONObject()
                    {{
                        put("groups", List.of(group));
                        put("value", normalizedExpression.get(group));
                    }});
                }
            }});
        }};
    }
}
