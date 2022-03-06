package util.Report;

import com.sun.tools.jconsole.JConsoleContext;
import com2pose.COM2POSE;
import org.json.JSONObject;
import org.json.JSONString;
import util.Configs.Configs;
import util.Configs.Modules.AbstractModule;
import util.FileManagement;
import util.MapSymbolAndEnsg;

import java.io.Console;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.lang.reflect.Field;
import java.net.URLEncoder;
import java.nio.charset.StandardCharsets;
import java.util.*;


public class PageGenerators
{
    static void generateGeneralPages() throws IOException
    {
        generateParameters();
        generateImportantLoci();
        generateTopLog2fc();
        generateCoOccurrence();
        generateOverview();
        generateRegressionPerformance();
    }

    static void generateImportantLoci() throws IOException
    {
        File sourceDir = COM2POSE.configs.igv.fileStructure.d_importantLoci.get();

        File targetDir = COM2POSE.configs.report.outputStructure.d_importantLoci.get();


        for (File groupFile : Objects.requireNonNull(sourceDir.listFiles()))
        {
            for (File imageFile : Objects.requireNonNull(groupFile.listFiles()))
            {
                if (!imageFile.getName().endsWith(".png"))
                {
                    continue;
                }

                for (Object importantTf : COM2POSE.configs.igv.importantLociAllPrioTf.get())
                {
                    if (imageFile.getName().startsWith(importantTf.toString()))
                    {
                        File targetFile = new File(
                                targetDir.getAbsolutePath() + File.separator + groupFile.getName() + File.separator +
                                        importantTf + File.separator + imageFile.getName());

                        FileManagement.copyFile(imageFile, targetFile);
                        break;
                    }
                }
            }
        }

        String important_loci = StructureElements.getFrame("Important loci",
                COM2POSE.configs.report.inputStructure.f_importantLoci.get());

        important_loci = important_loci.replace("{IMAGES}",
                StructureElements.generateImageSelector(targetDir.getName(), targetDir,
                        Arrays.asList(SelectorTypes.GROUPS, SelectorTypes.IMPORTANT_LOCI,
                                SelectorTypes.EMPTY_DROPDOWN)));

        FileManagement.writeHTML(COM2POSE.configs.report.outputStructure.f_importantLoci_html.get(), important_loci, 0);
    }

    static void generateTopLog2fc() throws IOException
    {
        File sourceDir = COM2POSE.configs.igv.fileStructure.d_igvTopLog2fc.get();

        File targetDir = COM2POSE.configs.report.outputStructure.d_topLog2fc.get();

        String top_log2fc =
                StructureElements.getFrame("Top log2fc", COM2POSE.configs.report.inputStructure.f_topLog2fc.get());

        top_log2fc = top_log2fc.replace("{TITLE}", "Top log2fc");

        for (File group_pairing : Objects.requireNonNull(sourceDir.listFiles()))
        {
            if (group_pairing.isDirectory())
            {
                for (File log2fc : Objects.requireNonNull(group_pairing.listFiles()))
                {
                    if (log2fc.isDirectory())
                    {
                        for (File imageFile : Objects.requireNonNull(log2fc.listFiles()))
                        {
                            if (imageFile.getName().endsWith(".png"))
                            {
                                File targetFile = new File(
                                        targetDir.getAbsolutePath() + File.separator + group_pairing.getName() +
                                                File.separator + log2fc.getName().split("_")[1] + File.separator +
                                                imageFile.getName());
                                FileManagement.copyFile(imageFile, targetFile);
                            }
                        }
                    }
                }
            }
        }

        top_log2fc = top_log2fc.replace("{IMAGES}",
                StructureElements.generateImageSelector(targetDir.getName(), targetDir,
                        Arrays.asList(SelectorTypes.GROUP_PAIRINGS, SelectorTypes.TOP_LOG2FC,
                                SelectorTypes.EMPTY_DROPDOWN), true, new JSONObject()));

        FileManagement.writeHTML(COM2POSE.configs.report.outputStructure.f_topLog2fc_html.get(), top_log2fc, 0);
    }

    static void generateCoOccurrence() throws IOException
    {
        File dataSource = COM2POSE.configs.distributionAnalysis.fileStructure.f_cooccurrence_frequencies.get();

        String input = FileManagement.loadFile(dataSource);
        String cooccurrence = StructureElements.getFrame("Co-Occurrence analysis",
                COM2POSE.configs.report.inputStructure.f_cooccurrence.get());

        Map<String, Map<String, Number>> data = new HashMap<>();

        ArrayList<String> tfs = new ArrayList<>();

        for (String entry : input.split("\n")[0].split("\t"))
        {
            if (!entry.isEmpty())
            {
                tfs.add(entry);
            }
        }

        for (String first : tfs)
        {
            data.put(first, new HashMap<>());
            for (String second : tfs)
            {
                String value = input.split("\n")[tfs.indexOf(first) + 1].split("\t")[tfs.indexOf(second) + 1];
                data.get(first).put(second, Double.parseDouble(value));
            }
        }

        cooccurrence = cooccurrence.replace("{TABLE}", StructureElements.getTabularData("coOccurrence", data));

        FileManagement.writeHTML(COM2POSE.configs.report.outputStructure.f_cooccurrence_html.get(), cooccurrence, 0);
    }

    static void generateOverview() throws IOException
    {
        File source = COM2POSE.configs.report.inputStructure.f_overview.get();
        File target = COM2POSE.configs.report.outputStructure.f_overview.get();

        String frame = StructureElements.getFrame("Overview", source);

        FileManagement.writeHTML(target, frame, 0);
    }

    static void generateParameters() throws IOException
    {
        Report.logger.logLine("[REPORT] Start generating report parameters page");
        String parameters =
                StructureElements.getFrame("Parameters", COM2POSE.configs.report.inputStructure.f_parameters.get());

        String moduleTemplate = FileManagement.loadFile(COM2POSE.configs.report.inputStructure.f_parameters_tool.get());

        StringBuilder sb_tools = new StringBuilder();

        JSONObject configs = COM2POSE.configs.getConfigsJSONObject(true);

        for (String moduleName : configs.keySet())
        {
            JSONObject module = configs.getJSONObject(moduleName);

            if (module.isEmpty())
            {
                continue;
            }

            String moduleString = moduleTemplate.replace("{TOOL_NAME}", moduleName);
            moduleString = moduleString.replace("{ID}", moduleName);

            StringBuilder sb_subModules = new StringBuilder();
            StringBuilder sb_simpleConfigs = new StringBuilder();

            for (String configName : module.keySet())
            {
                if (module.get(configName).getClass().equals(JSONObject.class))
                {
                    JSONObject subModule = module.getJSONObject(configName);

                    if (subModule.isEmpty())
                    {
                        continue;
                    }

                    String subModuleString = moduleTemplate.replace("{TOOL_NAME}", configName);
                    subModuleString = subModuleString.replace("{ID}", moduleName + "_" + configName);
                    subModuleString = subModuleString.replace("{SUBMODULES}", "");

                    StringBuilder sb_subConfigs = new StringBuilder();
                    for (String subConfigName : subModule.keySet())
                    {
                        sb_subConfigs.append("<p>");
                        sb_subConfigs.append(subConfigName);
                        sb_subConfigs.append(": ");
                        sb_subConfigs.append(subModule.get(subConfigName).toString());
                        sb_subConfigs.append("</p>");
                    }
                    subModuleString = subModuleString.replace("{SIMPLECONFIGS}", sb_subConfigs.toString());

                    sb_subModules.append(subModuleString);
                } else
                {
                    sb_simpleConfigs.append("<p>");
                    sb_simpleConfigs.append(configName);
                    sb_simpleConfigs.append(": ");
                    sb_simpleConfigs.append(module.get(configName).toString());
                    sb_simpleConfigs.append("</p>");
                }
            }

            moduleString = moduleString.replace("{SUBMODULES}", sb_subModules.toString());
            moduleString = moduleString.replace("{SIMPLECONFIGS}", sb_simpleConfigs.toString());

            sb_tools.append(moduleString);
        }

        parameters = parameters.replace("{TOOLS}", sb_tools.toString());

        //parameters = parameters.replace("{PARAMETERS}",
        //        FileManagement.loadFile(COM2POSE.configs.report.inputStructure.f_parameters_parameter.get()));

        FileManagement.writeHTML(COM2POSE.configs.report.outputStructure.f_parameters.get(), parameters, 0);
        Report.logger.logLine("[REPORT] Finished generating report parameters page");
    }

    static void generateHome(ArrayList<TranscriptionFactorGroup> transcriptionFactorGroups) throws IOException
    {
        Report.logger.logLine("[REPORT] Start generating report home page");
        String home = StructureElements.getFrame("Transcription factors",
                COM2POSE.configs.report.inputStructure.f_home.get());

        StringBuilder sb_tfs = new StringBuilder();


        int i = 1;

        for (TranscriptionFactorGroup tfGroup : transcriptionFactorGroups)
        {
            String prettyName = tfGroup.getName().replace("..", "_");
            sb_tfs.append("<div class='tf' id='" + prettyName + "'>");
            sb_tfs.append("<script>var targetGenes" + prettyName + " = " + tfGroup.getTargetGenes() + "</script>");

            if (!tfGroup.realGroup)
            {
                for (TranscriptionFactor transcriptionFactor : tfGroup.getTranscriptionFactors())
                {
                    String tf_string = FileManagement.loadFile(COM2POSE.configs.report.inputStructure.f_home_tf.get());

                    tf_string = tf_string.replace("{BUTTONBAR}", StructureElements.getButtonBar(tfGroup));

                    tf_string = tf_string.replace("{TF_NAME}", i + ". " + transcriptionFactor.getName());
                    tf_string = tf_string.replace("{ID}", String.valueOf(transcriptionFactor.getName().hashCode()));

                    tf_string = StructureElements.setBasicData(tf_string, transcriptionFactor);

                    tf_string = tf_string.replace("{GENEID}", transcriptionFactor.getGeneID());

                    tf_string = tf_string.replace("{HASVALIDATION}", tfGroup.hasValidation() ? "" : "disabled");
                    tf_string = tf_string.replace("{HASDISTRIBUTION}", tfGroup.hasDistribution() ? "" : "disabled");
                    tf_string = tf_string.replace("{HASREGRESSION}", tfGroup.hasRegression() ? "" : "disabled");

                    sb_tfs.append(tf_string);
                }
            } else
            {
                String tfGroupString =
                        FileManagement.loadFile(COM2POSE.configs.report.inputStructure.f_home_tfGroup.get());

                tfGroupString = tfGroupString.replace("{TF_NAME}", i + ". " + tfGroup.getName());
                tfGroupString = tfGroupString.replace("{ID}", String.valueOf(tfGroup.getName().hashCode()));
                tfGroupString = tfGroupString.replace("{BUTTONBAR}", StructureElements.getButtonBar(tfGroup));
                tfGroupString = tfGroupString.replace("{SINGLE_TFS}", StructureElements.getBasicData(tfGroup));

                sb_tfs.append(tfGroupString);
            }
            sb_tfs.append("</div>");

            i++;
        }
        home = home.replace("{TFS}", sb_tfs.toString());

        FileManagement.writeHTML(COM2POSE.configs.report.outputStructure.f_home.get(), home, 0);
        Report.logger.logLine("[REPORT] Finished generating report home page");
    }

    static void generateRegressionPerformance() throws IOException
    {
        {
            File parent = COM2POSE.configs.dynamite.fileStructure.d_output.get();

            for (File hm : Objects.requireNonNull(parent.listFiles()))
            {
                for (File group_pairing : Objects.requireNonNull(hm.listFiles()))
                {
                    FileManagement.copyFile(
                            new File(group_pairing.getAbsolutePath() + File.separator + "Performance_Barplots.png"),
                            new File(COM2POSE.configs.report.outputStructure.d_regression_performance_barplots.get()
                                    .getAbsolutePath() + File.separator + hm.getName() + File.separator +
                                    group_pairing.getName() + ".png"));

                    FileManagement.copyDirectory(group_pairing, new File(
                                    COM2POSE.configs.report.outputStructure.d_regression_performance_foldChanges.get()
                                            .getAbsolutePath() + File.separator + hm.getName() + File.separator +
                                            group_pairing.getName()), true,
                            "Misclassification_vs_Lambda_Fold_[1-9]_Integrated_Data_For_Classification.svg",
                            List.of("Misclassification_vs_Lambda_Fold", "Integrated_Data_For_Classification"));

                    FileManagement.copyFile(new File(group_pairing.getAbsolutePath() + File.separator +
                                    "Regression_Coefficients_Cross_Validation_Heatmap_Integrated_Data_For_Classification.svg"),
                            new File(COM2POSE.configs.report.outputStructure.d_regression_performance_heatmap.get()
                                    .getAbsolutePath() + File.separator + hm.getName() + File.separator +
                                    group_pairing.getName() + ".svg"));
                }
            }

        } // Copy files

        String frame = StructureElements.getFrame("Regression performance analysis",
                COM2POSE.configs.report.inputStructure.f_regressionPerformance.get());

        List<String> fileNames = new ArrayList<>();

        Report.existingValues.get(SelectorTypes.PERFORMANCE_CUTOFFS).forEach(cutoff -> fileNames.add(cutoff + ".svg"));

        File d_fold_changes = COM2POSE.configs.report.outputStructure.d_regression_performance_foldChanges.get();

        frame = frame.replace("{PERFORMANCE_FOLD_CHANGES}",
                StructureElements.generateImageSelector(d_fold_changes, "foldChanges",
                        Arrays.asList(Report.existingValues.get(SelectorTypes.HISTONE_MODIFICATIONS),
                                Report.existingValues.get(SelectorTypes.GROUP_PAIRINGS), fileNames), new JSONObject()));

        File d_barplots = COM2POSE.configs.report.outputStructure.d_regression_performance_barplots.get();
        frame = frame.replace("{PERFORMANCE_BARPLOT}", StructureElements.generateImageSelector("barplots", d_barplots,
                Arrays.asList(SelectorTypes.HISTONE_MODIFICATIONS, SelectorTypes.GROUP_PAIRINGS)));

        File d_heatmap = COM2POSE.configs.report.outputStructure.d_regression_performance_heatmap.get();
        fileNames.clear();
        Report.existingValues.get(SelectorTypes.GROUP_PAIRINGS)
                .forEach(group_pairing -> fileNames.add(group_pairing + ".svg"));

        frame = frame.replace("{PERFORMANCE_CROSS_VALIDATION_HEATMAP}",
                StructureElements.generateImageSelector(d_heatmap, "heatmap",
                        Arrays.asList(Report.existingValues.get(SelectorTypes.HISTONE_MODIFICATIONS), fileNames),
                        new JSONObject()));

        FileManagement.writeHTML(COM2POSE.configs.report.outputStructure.f_regression_performance_html.get(), frame, 2);
    }

    static boolean generateValidation(TranscriptionFactorGroup tfGroup) throws IOException
    {
        File templateFile = COM2POSE.configs.report.inputStructure.f_validation.get();

        File d_igv_screenshots = COM2POSE.configs.igv.fileStructure.d_root.get();

        File d_out_validation = COM2POSE.configs.report.outputStructure.d_validation.get();

        File d_heatmaps = COM2POSE.configs.distributionAnalysis.fileStructure.d_heatmaps.get();

        String frame = StructureElements.getFrame(tfGroup.getName() + " - Validation", templateFile);

        frame = StructureElements.setBasicData(frame, tfGroup);

        {
            File source = FileManagement.getFileIfInDirectory(d_heatmaps, tfGroup.getName(), false);

            String id = "validationHeatmaps";
            File target = new File(
                    d_out_validation.getAbsolutePath() + File.separator + tfGroup.getName() + File.separator + id);

            FileManagement.copyDirectory(source, target, false, ".+\\.png$", new ArrayList<>());

            JSONObject data = new JSONObject();

            {
                for (File d_hm : source.listFiles())
                {
                    JSONObject hmObject = new JSONObject();
                    for (File f_groupPairing : d_hm.listFiles())
                    {
                        if (f_groupPairing.getName().endsWith(".csv"))
                        {
                            String fileNameWithoutExtension = f_groupPairing.getName().replace(".csv", "");
                            hmObject.put(fileNameWithoutExtension, new JSONObject(new HashMap<>()
                            {{
                                put("content", FileManagement.loadFile(f_groupPairing));
                                put("fileExtension", ".csv");
                            }}));
                        }
                    }
                    data.put(d_hm.getName(), hmObject);
                }
            }

            frame = frame.replace("{VALIDATION_HEATMAP}", StructureElements.generateImageSelector(id, target,
                    Arrays.asList(SelectorTypes.HISTONE_MODIFICATIONS, SelectorTypes.GROUP_PAIRINGS), data));
        } // Heatmaps

        {
            File d_own_tf = COM2POSE.configs.igv.fileStructure.d_ownData.get();

            File source = FileManagement.getFileIfInDirectory(d_own_tf, "[0-9]+_" + tfGroup.getName(), false);

            String id = "validationOwnTF";
            File target = new File(
                    d_out_validation.getAbsolutePath() + File.separator + tfGroup.getName() + File.separator + id);

            source = null; // Deactivate generation of own tf chip-seq data

            if (source != null)
            {
                FileManagement.copyDirectory(source, target, true, ".+\\.png$", List.of(tfGroup.getName()));

                frame = frame.replace("{VALIDATION_OWN_TF}", StructureElements.generateImageSelector(id, target,
                        Arrays.asList(SelectorTypes.HISTONE_MODIFICATIONS, SelectorTypes.GROUPS,
                                SelectorTypes.EMPTY_DROPDOWN)));
            } else
            {
                frame = frame.replace("{VALIDATION_OWN_TF}", "");
            }

            frame = frame.replace("{VALIDATION_OWN_TF_DISABLED}", (source == null ? "disabled" : ""));
            frame = frame.replace("{VALIDATION_OWN_TF_VISIBLE}", (source == null ? "style='display: none'" : ""));
        } // Own tf

        {
            File d_chip_atlas = COM2POSE.configs.igv.fileStructure.d_chipAtlasData.get();

            File source = FileManagement.getFileIfInDirectory(d_chip_atlas, "[0-9]+_" + tfGroup.getName(), false);

            String id = "validationChipAtlas";

            File target = new File(
                    d_out_validation.getAbsolutePath() + File.separator + tfGroup.getName() + File.separator + id);

            if (source == null)
            {
                frame = frame.replace("{VALIDATION_CHIP_ATLAS}", "");
            } else
            {
                for (File d_group : Objects.requireNonNull(source.listFiles()))
                {
                    if (!d_group.isDirectory())
                    {
                        continue;
                    }
                    for (File d_hm : Objects.requireNonNull(d_group.listFiles()))
                    {
                        if (!d_hm.isDirectory())
                        {
                            continue;
                        }
                        for (File f_plot : Objects.requireNonNull(d_hm.listFiles()))
                        {
                            if (f_plot.getName().endsWith(".png"))
                            {
                                String relevantFileName = f_plot.getName();
                                relevantFileName = relevantFileName.replace("_" + tfGroup.getName(), "");

                                File targetFile = new File(
                                        target.getAbsolutePath() + File.separator + d_hm.getName() + File.separator +
                                                d_group.getName() + File.separator + relevantFileName);

                                FileManagement.copyFile(f_plot, targetFile);
                            }
                        }
                    }
                }
                String imageSelector = StructureElements.generateImageSelector(id, target,
                        Arrays.asList(SelectorTypes.HISTONE_MODIFICATIONS, SelectorTypes.GROUPS,
                                SelectorTypes.EMPTY_DROPDOWN));

                frame = frame.replace("{VALIDATION_CHIP_ATLAS}", imageSelector);
            }

            frame = frame.replace("{VALIDATION_CHIP_ATLAS_DISABLED}", target.exists() ? "" : "disabled");
            frame = frame.replace("{VALIDATION_CHIP_ATLAS_VISIBLE}", target.exists() ? "" : "style='display: none'");
        } // Chip Atlas

        {
            {
                File sourceDir = COM2POSE.configs.distributionAnalysis.fileStructure.d_logos_biophysicalModel.get();

                File targetFile = new File(
                        COM2POSE.configs.report.outputStructure.d_validation.get().getAbsolutePath() + File.separator +
                                tfGroup.getName() + File.separator +
                                COM2POSE.configs.report.outputStructure.s_validation_logos_biophysicalModel +
                                File.separator +
                                COM2POSE.configs.report.outputStructure.s_validation_logos_biophysicalModel_png);

                File tempDir = FileManagement.getFileIfInDirectory(sourceDir, "[0-9]+_" + tfGroup.getName(), false);

                File sourceFile = FileManagement.getFileIfInDirectory(tempDir, ".+\\.png$", true);
                File csvFile = FileManagement.getFileIfInDirectory(tempDir, ".+energy_matrix\\.csv$", true);

                if (sourceFile != null)
                {
                    FileManagement.copyFile(sourceFile, targetFile);

                    List<List<String>> options = new ArrayList<>();

                    options.add(new ArrayList<>(List.of(targetFile.getName())));

                    JSONObject data = new JSONObject();
                    data.put(targetFile.getName().replace(".png", ""), new JSONObject(new HashMap<>()
                    {{
                        put("content", FileManagement.loadFile(csvFile).replace("\t", ","));
                        put("fileExtension", ".csv");
                    }}));

                    frame = frame.replace("{BIOPHYSICAL_MODEL}",
                            StructureElements.generateImageSelector(targetFile.getParentFile(), "logosBiophysicalModel",
                                    options, data));
                }
            } // BIOPHYSICAL MODEL

            {
                boolean created = false;
                File sDir = COM2POSE.configs.distributionAnalysis.fileStructure.d_logos_tfSequence.get();

                File targetDir = new File(
                        COM2POSE.configs.report.outputStructure.d_validation.get().getAbsolutePath() + File.separator +
                                tfGroup.getName() + File.separator +
                                COM2POSE.configs.report.outputStructure.s_validation_logos_tfSequence.get());

                File sourceDir = FileManagement.getFileIfInDirectory(sDir, "[0-9]+_" + tfGroup.getName(), false);

                if (sourceDir != null)
                {
                    ArrayList<String> files = new ArrayList<>();
                    for (File file : Objects.requireNonNull(sourceDir.listFiles()))
                    {
                        if (file.getName().endsWith(".svg"))
                        {
                            File targetFile = new File(targetDir.getAbsolutePath() + File.separator + file.getName());
                            FileManagement.copyFile(file, targetFile);
                            files.add(targetFile.getName());
                        }
                    }

                    if (files.size() > 0)
                    {
                        List<List<String>> options = new ArrayList<>();
                        options.add(files);

                        JSONObject data = new JSONObject();

                        for (String imageFileName : files)
                        {
                            String nameWithoutExtension = imageFileName.replace(".svg", "");

                            File jsonFile =
                                    FileManagement.getFileIfInDirectory(sourceDir, nameWithoutExtension + ".json",
                                            true);

                            if (jsonFile != null)
                            {
                                data.put(nameWithoutExtension, new JSONObject(new HashMap<>()
                                {{
                                    put("content", FileManagement.loadFile(jsonFile));
                                    put("fileExtension", ".json");
                                }}));
                            }
                        }


                        frame = frame.replace("{TF_SEQUENCE}",
                                StructureElements.generateImageSelector(targetDir, "logosTfSequence", options, data));
                        created = true;
                    }
                }
                if (!created)
                {
                    frame = frame.replace("{TF_SEQUENCE}", "");
                }
                frame = frame.replace("{TF_SEQUENCE_DISABLED}", created ? "" : "disabled");
            } // TF Sequence

            {
                File sDir = COM2POSE.configs.distributionAnalysis.fileStructure.d_logos_tfBindingSequence.get();

                File targetDir = new File(
                        COM2POSE.configs.report.outputStructure.d_validation.get().getAbsolutePath() + File.separator +
                                tfGroup.getName() + File.separator +
                                COM2POSE.configs.report.outputStructure.s_validation_logos_tfBindingSequence.get());

                File sourceDir = FileManagement.getFileIfInDirectory(sDir, "[0-9]+_" + tfGroup.getName(), false);

                if (sourceDir != null)
                {
                    FileManagement.copyDirectory(sourceDir, targetDir, true, ".*_biophysical_model.png$",
                            Arrays.asList(tfGroup.getName(), "biophysical_model"));

                    frame = frame.replace("{TF_BINDING_SEQUENCE}",
                            StructureElements.generateImageSelector("logosTfBindingSequence", targetDir,
                                    List.of(SelectorTypes.DISTRIBUTION_OPTIONS)));
                }

                frame = frame.replace("{TF_BINDING_SEQUENCE_VISIBLE}",
                        sourceDir == null ? "style='display: none'" : "");

                if (sourceDir != null)
                {

                }
            } // TF binding sequence
        } // LOGOS

        {
            File sourceDir = COM2POSE.configs.igv.fileStructure.d_igvDcgTargetGenes.get();

            File targetDir = COM2POSE.configs.report.outputStructure.d_validation.get();

            File source = FileManagement.getFileIfInDirectory(sourceDir, tfGroup.getName(), false);

            String id = "validationIGV";
            File target =
                    new File(targetDir.getAbsolutePath() + File.separator + tfGroup.getName() + File.separator + id);

            frame = frame.replace("{VALIDATION_IGV_DISABLED}", (source == null) ? "disabled" : "");

            if (source != null)
            {
                FileManagement.copyDirectory(source, target, true, ".+\\.png$", new ArrayList<>());
            }

            JSONObject data = new JSONObject();
            {
                Map<String, Set<String>> csvSetData = new HashMap<>();

                for (File d_hm : target.listFiles())
                {
                    csvSetData.put(d_hm.getName(), new HashSet<>());
                    for (File d_groupPairing : d_hm.listFiles())
                    {
                        for (File f_targetGene : d_groupPairing.listFiles())
                        {
                            String name = f_targetGene.getName();
                            if (name.endsWith(".png"))
                            {
                                name = name.split("_")[1];
                                name = name.split("\\.")[0];
                                csvSetData.get(d_hm.getName()).add(name);
                            }
                        }
                    }
                }

                for (String hm : Report.existingValues.get(SelectorTypes.HISTONE_MODIFICATIONS))
                {
                    StringBuilder sb_hm = new StringBuilder();
                    sb_hm.append("GeneSymbol,ENSG\n");

                    for (String targetGeneSymbol : csvSetData.get(hm))
                    {
                        sb_hm.append(targetGeneSymbol);
                        sb_hm.append(",");
                        try
                        {
                            sb_hm.append(MapSymbolAndEnsg.symbolToEnsg(targetGeneSymbol));
                        } catch (NoSuchFieldException | FileNotFoundException e)
                        {
                            Report.logger.warn(e.getMessage());
                        }
                        sb_hm.append("\n");
                    }

                    data.put(hm, new JSONObject(new HashMap<>()
                    {{
                        put("content", URLEncoder.encode(sb_hm.toString(), StandardCharsets.UTF_8));
                        put("fileExtension", ".csv");
                    }}));
                }
            }

            frame = frame.replace("{VALIDATION_IGV}", (source == null) ? "" :
                    StructureElements.generateImageSelector(id, target,
                            Arrays.asList(SelectorTypes.HISTONE_MODIFICATIONS, SelectorTypes.GROUP_PAIRINGS,
                                    SelectorTypes.EMPTY_DROPDOWN), true, data));
        } // IGV

        frame = frame.replace("{TFNAME}", tfGroup.getName());

        frame = StructureElements.setGeneCardLinks(frame, tfGroup);

        FileManagement.writeHTML(new File(
                COM2POSE.configs.report.outputStructure.d_validation.get().getAbsolutePath() + File.separator +
                        tfGroup.getName() + File.separator + tfGroup.getName() + ".html"), frame, 2);

        return true;
    }

    static boolean generateDistribution(TranscriptionFactorGroup tfGroup) throws IOException
    {
        String id = "distributionPlots";

        File templateFile = COM2POSE.configs.report.inputStructure.f_distribution.get();

        File d_distribution_output = new File(
                COM2POSE.configs.report.outputStructure.d_distribution.get().getAbsolutePath() + File.separator +
                        tfGroup.getName());

        File d_distribution_plots = new File(d_distribution_output.getAbsolutePath() + File.separator + id);


        File d_plots_hm = COM2POSE.configs.distributionAnalysis.fileStructure.d_plots_hm.get();

        for (File d_hm : Objects.requireNonNull(d_plots_hm.listFiles()))
        {
            for (File f_plot : Objects.requireNonNull(d_hm.listFiles()))
            {
                if (f_plot.getName().substring(0, f_plot.getName().lastIndexOf(".")).equals(tfGroup.getName()))
                {
                    FileManagement.copyFile(f_plot, new File(
                            d_distribution_plots.getAbsolutePath() + File.separator + d_hm.getName() + ".png"));
                    break;
                }
            }
        }

        String frame = StructureElements.getFrame(tfGroup.getName() + " - Distribution", templateFile);

        frame = frame.replace("{TFNAME}", tfGroup.getName());
        frame = StructureElements.setBasicData(frame, tfGroup);


        frame = frame.replace("{DISTRIBUTION_PLOTS}", StructureElements.generateImageSelector(id, d_distribution_plots,
                List.of(SelectorTypes.DISTRIBUTION_OPTIONS)));


        frame = StructureElements.setGeneCardLinks(frame, tfGroup);


        {
            File hmsDir = COM2POSE.configs.distributionAnalysis.fileStructure.d_stats_hm.get();
            HashMap<String, Integer> ranks = new HashMap<>();
            HashMap<String, Integer> sizes = new HashMap<>();

            for (File hmDir : Objects.requireNonNull(hmsDir.listFiles()))
            {
                if (!hmDir.isDirectory())
                {
                    continue;
                }

                File statsFile = new File(hmDir.getAbsolutePath() + File.separator +
                        COM2POSE.configs.distributionAnalysis.fileStructure.s_stats_csv);

                if (statsFile.exists())
                {
                    try
                    {
                        ranks.put(hmDir.getName(), Integer.parseInt(
                                FileManagement.findValueInTable(tfGroup.getName(), 1, 0, statsFile, "\t", true)));
                    } catch (NoSuchFieldException e)
                    {
                        ranks.put(hmDir.getName(), -1);
                    } finally
                    {
                        String[] lines = FileManagement.loadFile(statsFile).split("\n");
                        String lastLine = lines[lines.length - 1];
                        sizes.put(hmDir.getName(), Integer.parseInt(lastLine.split("\t")[0]));
                    }
                }
            }

            StringBuilder sb_dcg = new StringBuilder();
            sb_dcg.append("<table>");
            sb_dcg.append("<tr>");
            sb_dcg.append("<th>Modification</th>");
            sb_dcg.append("<th>Availability</th>");
            sb_dcg.append("<th>Rank</th>");
            sb_dcg.append("</tr>");

            for (Map.Entry<String, Integer> entry : ranks.entrySet())
            {
                sb_dcg.append("<tr>");
                sb_dcg.append("<td>");
                sb_dcg.append(entry.getKey());
                sb_dcg.append("</td>");
                sb_dcg.append("<td>");
                sb_dcg.append("<image style='height: 30px' src='");
                sb_dcg.append(entry.getValue() == -1 ? "{RELATIVATION}MEDIA/not_available.png" :
                        "{RELATIVATION" + "}MEDIA/is_available.png");
                sb_dcg.append("'>");
                sb_dcg.append("</td>");
                sb_dcg.append("<td>");
                sb_dcg.append(entry.getValue() > 0 ? entry.getValue() : "-");
                sb_dcg.append("/");
                sb_dcg.append(sizes.get(entry.getKey()));
                sb_dcg.append("</td>");
                sb_dcg.append("</tr>");
            }
            sb_dcg.append("</table>");

            Double score;
            try
            {
                score = Double.parseDouble(FileManagement.findValueInTable(tfGroup.getName(), 1, 2,
                        COM2POSE.configs.distributionAnalysis.fileStructure.f_dcg_stats.get(), "\t", true));
                frame = frame.replace("{DCG_SCORE}", Report.formatter.format(score));
            } catch (NoSuchFieldException ignored)
            {
                frame = frame.replace("{DCG_SCORE}", "-");
            }

            frame = frame.replace("{DISCOUNTED_CUMULATIVE_GAIN}", sb_dcg.toString());
        }

        FileManagement.writeHTML(new File(d_distribution_output + File.separator + tfGroup.getName() + ".html"), frame,
                2);

        return true;
    }

    static boolean generateRegression(TranscriptionFactorGroup tfGroup) throws IOException
    {
        File templateFile = COM2POSE.configs.report.inputStructure.f_regression.get();

        File d_in_plots = COM2POSE.configs.plots.fileStructure.d_output.get();

        File d_out_regression = COM2POSE.configs.report.outputStructure.d_regression.get();

        String frame = StructureElements.getFrame(tfGroup.getName() + " - Regression", templateFile);

        frame = frame.replace("{TFNAME}", tfGroup.getName());

        File tfDir = new File(d_out_regression.getAbsolutePath() + File.separator + tfGroup.getName());

        String overviewHeatmapsID = "overviewHeatmaps";
        String overviewCoefficientsID = "overviewCoefficients";

        for (File hm_dir : Objects.requireNonNull(d_in_plots.listFiles()))
        {
            if (!hm_dir.isDirectory())
            {
                continue;
            }

            for (File cutoff_dir : Objects.requireNonNull((hm_dir.listFiles())))
            {
                if (!cutoff_dir.isDirectory())
                {
                    continue;
                }

                for (File image_file : Objects.requireNonNull(cutoff_dir.listFiles()))
                {
                    if (!image_file.getName().endsWith(".png"))
                    {
                        continue;
                    }

                    String relevantFileName = image_file.getName();
                    List<String> removables =
                            Arrays.asList(".png", hm_dir.getName(), "threshold", cutoff_dir.getName());
                    for (String removable : removables)
                    {
                        relevantFileName = relevantFileName.replace(removable, "");
                    }
                    while (relevantFileName.endsWith("_"))
                    {
                        relevantFileName = relevantFileName.substring(0, relevantFileName.length() - 1);
                    }

                    while (relevantFileName.startsWith("_"))
                    {
                        relevantFileName = relevantFileName.substring(1);
                    }

                    String targetCategory =
                            relevantFileName.contains("stage") ? overviewHeatmapsID : overviewCoefficientsID;

                    File targetFile = new File(
                            tfDir.getAbsolutePath() + File.separator + targetCategory + File.separator +
                                    hm_dir.getName() + File.separator + cutoff_dir.getName() + File.separator +
                                    relevantFileName + ".png");

                    FileManagement.copyFile(image_file, targetFile);
                }
            }
        }

        String overviewCoefficients = StructureElements.generateImageSelector(overviewCoefficientsID,
                new File(tfDir.getAbsolutePath() + File.separator + overviewCoefficientsID),
                Arrays.asList(SelectorTypes.HISTONE_MODIFICATIONS, SelectorTypes.REGRESSION_CUTOFFS,
                        SelectorTypes.EMPTY_DROPDOWN));

        String overviewHeatmaps = StructureElements.generateImageSelector(overviewHeatmapsID,
                new File(tfDir.getAbsolutePath() + File.separator + overviewHeatmapsID),
                Arrays.asList(SelectorTypes.HISTONE_MODIFICATIONS, SelectorTypes.REGRESSION_CUTOFFS,
                        SelectorTypes.EMPTY_DROPDOWN));

        frame = frame.replace("{COEFFICIENTS}",
                StructureElements.getTabularData("regressionCoefficients", tfGroup.regressionCoefficients));

        frame = frame.replace("{OVERVIEW_COEFFICIENTS}", overviewCoefficients);

        frame = frame.replace("{OVERVIEW_HEATMAPS}", overviewHeatmaps);
        frame = StructureElements.setBasicData(frame, tfGroup);

        frame = StructureElements.setGeneCardLinks(frame, tfGroup);

        FileManagement.writeHTML(new File(
                COM2POSE.configs.report.outputStructure.d_regression.get().getAbsolutePath() + File.separator +
                        tfGroup.getName() + File.separator + tfGroup.getName() + ".html"), frame, 2);

        return true;
    }
}
