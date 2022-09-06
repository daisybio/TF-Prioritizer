package lib.AngularReport;

import lib.ExecutableStep;
import org.json.JSONObject;
import tfprio.tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.FileFilters.Filters;
import util.Logger;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ExecutorService;

import static util.FileManagement.*;
import static util.ScriptExecution.executeAndWait;

public class Generate extends ExecutableStep
{
    private final AbstractConfig<File> f_input_transcriptionFactors =
            TFPRIO.configs.distributionAnalysis.fileStructure.f_dcg_stats;
    private final AbstractConfig<File> d_input_confusionMatrixes =
            TFPRIO.configs.confusionMatrixes.fileStructure.d_output;
    private final AbstractConfig<File> d_input_report = TFPRIO.configs.angularReport.d_angularSource;
    private final GeneratedFileStructure f_output_reportJson = TFPRIO.configs.angularReport.fileStructure.f_data;
    private final GeneratedFileStructure d_output_data = TFPRIO.configs.angularReport.fileStructure.d_data;
    private final GeneratedFileStructure d_preprocessing = TFPRIO.configs.angularReport.fileStructure.d_preprocessing;
    private final GeneratedFileStructure d_output = TFPRIO.configs.angularReport.fileStructure.d_output;

    JSONObject importantLoci;

    JSONObject topLog2fc;
    JSONObject statisticalEvaluation;

    ArrayList<JSONObject> coOccurrence;

    public static void linkFiles(JSONObject map, File d_target, ExecutorService executorService, Logger logger)
    {
        Set<String> keysToRemove = new HashSet<>();

        for (String key : map.keySet())
        {
            File newTarget = extend(d_target, key);

            Object value = map.get(key);

            if (value.getClass().equals(String.class))
            {
                File sourceFile = new File((String) value);
                String extension = sourceFile.getName().substring(sourceFile.getName().lastIndexOf("."));

                File targetFile = new File(newTarget.getAbsolutePath() + extension);

                executorService.submit(() -> softLink(targetFile, sourceFile, logger));

                map.put(key, targetFile.getAbsolutePath()
                        .replace(TFPRIO.configs.angularReport.fileStructure.d_data.get().getAbsolutePath(), ""));

            } else if (value.getClass().equals(JSONObject.class))
            {
                linkFiles((JSONObject) value, newTarget, executorService, logger);
            } else if (value.getClass().equals(java.math.BigDecimal.class) ||
                    value.getClass().equals(java.lang.Integer.class))
            {
                continue;
            } else
            {
                logger.error("Illegal value type detected: " + value.getClass());
            }
        }
    }

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(List.of(f_input_transcriptionFactors));
    }

    @Override public Set<GeneratedFileStructure> getCreatedFileStructure()
    {
        return new HashSet<>(Arrays.asList(f_output_reportJson, d_output_data));
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>();
    }

    @Override protected void execute()
    {
        copyDirectory(d_input_report.get(), d_preprocessing.get(), false);

        collectImportantLoci();
        collectTopLog2fc();
        collectCoOccurrence();
        collectStatisticalEvaluation();
        logger.info("Finished collecting general data.");

        List<TranscriptionFactorGroup> transcriptionFactorGroups = new ArrayList<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(f_input_transcriptionFactors.get())))
        {
            String inputLine;
            reader.readLine();

            while ((inputLine = reader.readLine()) != null)
            {
                String[] split = inputLine.split("\t");

                if (split.length != 3)
                {
                    logger.warn("Irregular line length in tf ranking file detected");
                }

                String symbol = split[1];

                TranscriptionFactorGroup transcriptionFactorGroup =
                        new TranscriptionFactorGroup(symbol, logger, executorService);
                transcriptionFactorGroups.add(transcriptionFactorGroup);
                transcriptionFactorGroup.collectData();
            }
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }
        finishAllQueuedThreads();

        saveJson(transcriptionFactorGroups);

        String script = "cd " + d_preprocessing.get().getAbsolutePath() +
                "\nrm -rf node_modules ; rm -f package-lock.json ; npm install --include-dev\nng build --output-path=" +
                d_output.get().getAbsolutePath();

        try
        {
            writeFile(TFPRIO.configs.angularReport.fileStructure.f_script.get(), script);
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        executeAndWait(TFPRIO.configs.angularReport.fileStructure.f_script.get(), logger);
    }

    private void collectImportantLoci()
    {
        File d_source = TFPRIO.configs.igv.fileStructure.d_importantLoci.get();

        Map<String, Map<String, Map<String, Map<String, File>>>> group_importantLocus_targetGene_fileType_file =
                new HashMap<>();

        for (String group : TFPRIO.groupsToHms.keySet())
        {
            File d_hm = extend(d_source, group);

            group_importantLocus_targetGene_fileType_file.put(group, new HashMap<>());

            for (File f_plot : Objects.requireNonNull(d_hm.listFiles(Filters.getSuffixFilter(".png"))))
            {
                String targetGene = f_plot.getName().substring(0, f_plot.getName().lastIndexOf("."));

                for (String importantLocus : TFPRIO.configs.igv.importantLociAllPrioTf.get())
                {
                    if (targetGene.contains(importantLocus))
                    {
                        if (!group_importantLocus_targetGene_fileType_file.get(group).containsKey(importantLocus))
                        {
                            group_importantLocus_targetGene_fileType_file.get(group)
                                    .put(importantLocus, new HashMap<>());
                        }

                        group_importantLocus_targetGene_fileType_file.get(group).get(importantLocus)
                                .put(targetGene, new HashMap<>());

                        group_importantLocus_targetGene_fileType_file.get(group).get(importantLocus).get(targetGene)
                                .put("plot", f_plot);
                    }
                }
            }
        }

        importantLoci = new JSONObject(group_importantLocus_targetGene_fileType_file);

        linkFiles(importantLoci, extend(d_output_data.get(), "importantLoci"), executorService, logger);
    }

    private void collectTopLog2fc()
    {
        File d_source = TFPRIO.configs.igv.fileStructure.d_igvTopLog2fc.get();

        Map<String, Map<String, Map<String, Map<String, File>>>> groupPairing_regulationType_targetGene_fileType_file =
                new HashMap<>();


        for (String groupPairing : TFPRIO.groupCombinationsToHms.keySet())
        {
            groupPairing_regulationType_targetGene_fileType_file.put(groupPairing, new HashMap<>());

            File d_groupPairing = extend(d_source, groupPairing);

            if (!d_groupPairing.exists())
            {
                continue;
            }

            for (File d_regulationType : Objects.requireNonNull(d_groupPairing.listFiles(Filters.directoryFilter)))
            {
                String regulationType = d_regulationType.getName().replaceAll("\\d+_", "");

                groupPairing_regulationType_targetGene_fileType_file.get(groupPairing)
                        .put(regulationType, new HashMap<>());

                for (File f_targetGenePlot : Objects.requireNonNull(
                        d_regulationType.listFiles(Filters.getSuffixFilter(".png"))))
                {
                    String targetGene = f_targetGenePlot.getName().replaceAll("\\d+_", "");
                    targetGene = targetGene.substring(0, targetGene.lastIndexOf("."));

                    groupPairing_regulationType_targetGene_fileType_file.get(groupPairing).get(regulationType)
                            .put(targetGene, new HashMap<>()
                            {{
                                put("plot", f_targetGenePlot);
                            }});
                }
            }
        }
        topLog2fc = new JSONObject(groupPairing_regulationType_targetGene_fileType_file);

        linkFiles(topLog2fc, extend(d_output_data.get(), "topLog2fc"), executorService, logger);
    }

    private void collectCoOccurrence()
    {
        File f_source = TFPRIO.configs.distributionAnalysis.fileStructure.f_cooccurrence_frequencies.get();

        List<String> transcriptionFactors = new ArrayList<>();

        Map<String, Map<String, Double>> tf_tf_value = new HashMap<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(f_source)))
        {
            String firstLine = reader.readLine();

            for (String transcriptionFactor : firstLine.split("\t"))
            {
                if (!transcriptionFactor.isBlank())
                {
                    transcriptionFactors.add(transcriptionFactor);
                }
            }

            String inputLine;
            while ((inputLine = reader.readLine()) != null)
            {
                String[] split = inputLine.split("\t");
                String rowTfName = split[0];

                Map<String, Double> tf_value = new HashMap<>();

                for (int i = 1; i < split.length; i++)
                {
                    double value = Double.parseDouble(split[i]);
                    String columnTfName = transcriptionFactors.get(i - 1);

                    tf_value.put(columnTfName, value);
                }

                tf_tf_value.put(rowTfName, tf_value);
            }
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        coOccurrence = new ArrayList<>()
        {{
            for (String tf1 : tf_tf_value.keySet())
            {
                for (String tf2 : tf_tf_value.get(tf1).keySet())
                {
                    if (tf1.compareTo(tf2) >= 0)
                    {
                        continue;
                    }

                    add(new JSONObject()
                    {{
                        put("groups", Arrays.asList(tf1, tf2));
                        put("value", tf_tf_value.get(tf1).get(tf2));
                    }});
                }
            }
        }};
    }

    private void collectStatisticalEvaluation()
    {
        statisticalEvaluation = new JSONObject();

        for (File f_tf : Objects.requireNonNull(
                d_input_confusionMatrixes.get().listFiles(Filters.getSuffixFilter(".json"))))
        {
            String tfName = f_tf.getName().replace(".json", "");
            String content = readFile(f_tf, logger);
            assert content != null;
            JSONObject tfJsonObject = new JSONObject(content);
            statisticalEvaluation.put(tfName, tfJsonObject);
        }
        statisticalEvaluation.put("metricsPlot",
                extend(d_input_confusionMatrixes.get(), "metrics.png").getAbsolutePath());
        linkFiles(statisticalEvaluation, extend(d_output_data.get(), "metrics"), executorService, logger);
    }

    private void saveJson(List<TranscriptionFactorGroup> transcriptionFactorGroups)
    {
        JSONObject jsonObject = new JSONObject()
        {{
            put("configs", TFPRIO.configs.getConfigsJSONObject(true, true));
            put("importantLoci", importantLoci);
            put("topLog2fc", topLog2fc);
            put("coOccurrenceAnalysis", coOccurrence);
            put("backgroundStats", TranscriptionFactorGroup.getBackgroundStats());
            put("statisticalEvaluation", statisticalEvaluation);
            put("existingValues", new JSONObject(new HashMap<>()
            {{
                put("hm", TFPRIO.existingHms);
                put("group", TFPRIO.groupsToHms.keySet());
                put("groupPairing", TFPRIO.groupCombinationsToHms.keySet());
            }}));
            put("transcriptionFactorGroups", new ArrayList<>()
            {{
                transcriptionFactorGroups.forEach(
                        transcriptionFactorGroup -> add(transcriptionFactorGroup.toJSONObject()));
            }});
        }};


        try
        {
            writeFile(f_output_reportJson.get(), "export const tfData=" + jsonObject.toString(4));
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }
    }
}