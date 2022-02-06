package util.Report;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.TimeUnit;

public class PageGenerators
{
    static void generateGeneralPages() throws IOException
    {
        generateParameters();
        generateImportantLoci();
        generateTopLog2fc();
        generateCoOccurrence();
        generateOverview();
    }

    static void generateImportantLoci() throws IOException
    {
        File sourceDir = new File(Report.options_intern.com2pose_working_directory + File.separator +
                Report.options_intern.folder_out_igv + File.separator +
                Report.options_intern.folder_out_igv_important_loci);

        File targetDir = new File(Report.options_intern.com2pose_working_directory + File.separator +
                Report.options_intern.d_out_important_loci);


        for (File groupFile : Objects.requireNonNull(sourceDir.listFiles()))
        {
            for (File imageFile : Objects.requireNonNull(groupFile.listFiles()))
            {
                if (!imageFile.getName().endsWith(".png"))
                {
                    continue;
                }

                for (String importantTf : Report.options_intern.igv_important_locus_all_prio_tf)
                {
                    if (imageFile.getName().startsWith(importantTf))
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
                Report.options_intern.path_to_COM2POSE + File.separator +
                        Report.options_intern.f_report_resources_important_loci_html);

        important_loci = important_loci.replace("{IMAGES}",
                StructureElements.generateImageSelector(targetDir.getName(), targetDir,
                        Arrays.asList(SelectorTypes.GROUPS, SelectorTypes.IMPORTANT_LOCI,
                                SelectorTypes.EMPTY_DROPDOWN)));

        FileManagement.writeHTML(Report.options_intern.com2pose_working_directory + File.separator +
                Report.options_intern.f_out_report_important_loci_html, important_loci, 0);
    }

    static void generateTopLog2fc() throws IOException
    {
        File sourceDir = new File(Report.options_intern.com2pose_working_directory + File.separator +
                Report.options_intern.folder_out_igv + File.separator +
                Report.options_intern.folder_out_igv_top_log2fc);

        File targetDir = new File(Report.options_intern.com2pose_working_directory + File.separator +
                Report.options_intern.d_out_top_log2fc);

        String top_log2fc = StructureElements.getFrame("Top log2fc",
                Report.options_intern.path_to_COM2POSE + File.separator +
                        Report.options_intern.f_report_resources_top_log2fc_html);

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
                                SelectorTypes.EMPTY_DROPDOWN)));

        FileManagement.writeHTML(Report.options_intern.com2pose_working_directory + File.separator +
                Report.options_intern.f_out_report_top_log2fc_html, top_log2fc, 0);
    }

    static void generateCoOccurrence() throws IOException
    {
        File dataSource = new File(Report.options_intern.com2pose_working_directory + File.separator +
                Report.options_intern.folder_out_distribution + File.separator +
                Report.options_intern.folder_out_distribution_cooccurence + File.separator +
                Report.options_intern.file_suffix_cooccurence_frequencies);

        String input = FileManagement.loadFile(dataSource.getAbsolutePath());
        String cooccurrence = StructureElements.getFrame("Co-Occurrence analysis",
                Report.options_intern.path_to_COM2POSE + File.separator +
                        Report.options_intern.f_report_resources_cooccurrence_html);

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

        FileManagement.writeHTML(Report.options_intern.com2pose_working_directory + File.separator +
                Report.options_intern.f_out_report_cooccurrence_html, cooccurrence, 0);
    }

    static void generateOverview() throws IOException
    {
        File source = new File(Report.options_intern.path_to_COM2POSE + File.separator +
                Report.options_intern.f_report_resources_overview_html);
        File target = new File(Report.options_intern.com2pose_working_directory + File.separator +
                Report.options_intern.f_out_report_overview);

        String frame = StructureElements.getFrame("Overview", source.getAbsolutePath());

        FileManagement.writeHTML(target.getAbsolutePath(), frame, 0);
    }

    static void generateParameters() throws IOException
    {
        Report.logger.logLine("[REPORT] Start generating report parameters page");
        String parameters = StructureElements.getFrame("Parameters",
                Report.options_intern.path_to_COM2POSE + File.separator +
                        Report.options_intern.f_report_resources_parameters_parameters_html);

        parameters = parameters.replace("{TOOLS}", FileManagement.loadFile(
                Report.options_intern.path_to_COM2POSE + File.separator +
                        Report.options_intern.f_report_resources_parameters_tool_html));
        parameters = parameters.replace("{PARAMETERS}", FileManagement.loadFile(
                Report.options_intern.path_to_COM2POSE + File.separator +
                        Report.options_intern.f_report_resources_parameters_parameter_html));

        FileManagement.writeHTML(Report.options_intern.com2pose_working_directory + File.separator +
                Report.options_intern.f_out_report_parameters, parameters, 0);
        Report.logger.logLine("[REPORT] Finished generating report parameters page");
    }

    static void generateHome(ArrayList<TranscriptionFactorGroup> transcriptionFactorGroups) throws IOException
    {
        Report.logger.logLine("[REPORT] Start generating report home page");
        String home = StructureElements.getFrame("Transcription factors",
                Report.options_intern.path_to_COM2POSE + File.separator +
                        Report.options_intern.f_report_resources_home_home_html);

        StringBuilder sb_tfs = new StringBuilder();


        int i = 1;

        for (TranscriptionFactorGroup tfGroup : transcriptionFactorGroups)
        {
            if (!tfGroup.realGroup)
            {
                for (TranscriptionFactor transcriptionFactor : tfGroup.getTranscriptionFactors())
                {
                    String tf_string = FileManagement.loadFile(Report.options_intern.path_to_COM2POSE + File.separator +
                            Report.options_intern.f_report_resources_home_tf_html);

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
                String tfGroupString = FileManagement.loadFile(Report.options_intern.path_to_COM2POSE + File.separator +
                        Report.options_intern.f_report_resources_home_tfGroup_html);

                tfGroupString = tfGroupString.replace("{TF_NAME}", i + ". " + tfGroup.getName());
                tfGroupString = tfGroupString.replace("{ID}", String.valueOf(tfGroup.getName().hashCode()));
                tfGroupString = tfGroupString.replace("{BUTTONBAR}", StructureElements.getButtonBar(tfGroup));
                tfGroupString = tfGroupString.replace("{SINGLE_TFS}", StructureElements.getBasicData(tfGroup));

                sb_tfs.append(tfGroupString);
            }

            i++;
        }
        home = home.replace("{TFS}", sb_tfs.toString());

        FileManagement.writeHTML(Report.options_intern.com2pose_working_directory + File.separator +
                Report.options_intern.f_out_report_home, home, 0);
        Report.logger.logLine("[REPORT] Finished generating report home page");
    }

    static boolean generateValidation(TranscriptionFactorGroup tfGroup) throws IOException
    {
        File templateFile = new File(Report.options_intern.path_to_COM2POSE + File.separator +
                Report.options_intern.f_report_resources_validation_validation_html);

        File d_igv_screenshots = new File(Report.options_intern.com2pose_working_directory + File.separator +
                Report.options_intern.folder_out_igv);

        File d_out_validation = new File(Report.options_intern.com2pose_working_directory + File.separator +
                Report.options_intern.d_out_validation);

        File d_heatmaps = new File(Report.options_intern.com2pose_working_directory + File.separator +
                Report.options_intern.folder_out_heatmap);

        String frame = StructureElements.getFrame(tfGroup.getName() + " - Validation", templateFile.getAbsolutePath());

        frame = frame.replace("{TFNAME}", tfGroup.getName());

        frame = StructureElements.setBasicData(frame, tfGroup);

        {
            File source = FileManagement.getFileIfInDirectory(d_heatmaps, tfGroup.getName(), false);

            String id = "validationHeatmaps";
            File target = new File(
                    d_out_validation.getAbsolutePath() + File.separator + tfGroup.getName() + File.separator + id);

            FileManagement.copyDirectory(source, target, false, ".+\\.png$", new ArrayList<>());

            frame = frame.replace("{VALIDATION_HEATMAP}", StructureElements.generateImageSelector(id, target,
                    Arrays.asList(SelectorTypes.HISTONE_MODIFICATIONS, SelectorTypes.GROUP_PAIRINGS)));
        } // Heatmaps

        {
            File d_own_tf = new File(d_igv_screenshots.getAbsolutePath() + File.separator +
                    Report.options_intern.folder_out_igv_own_data);

            File source = FileManagement.getFileIfInDirectory(d_own_tf, "[0-9]+_" + tfGroup.getName(), false);

            String id = "validationOwnTF";
            File target = new File(
                    d_out_validation.getAbsolutePath() + File.separator + tfGroup.getName() + File.separator + id);

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
        } // Own tf

        {
            File d_chip_atlas = new File(d_igv_screenshots.getAbsolutePath() + File.separator +
                    Report.options_intern.folder_out_igv_chip_atlas_data);

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
        } // Chip Atlas

        {
            File sourceParentDir = new File(Report.options_intern.com2pose_working_directory + File.separator +
                    Report.options_intern.folder_out_distribution + File.separator +
                    Report.options_intern.folder_out_distribution_logos);

            {
                File sourceDir = new File(sourceParentDir.getAbsolutePath() + File.separator +
                        Report.options_intern.folder_out_distribution_logos_biophysiscal_model);

                File targetFile = new File(Report.options_intern.com2pose_working_directory + File.separator +
                        Report.options_intern.d_out_validation + File.separator + tfGroup.getName() + File.separator +
                        Report.options_intern.d_out_validation_logos_biophysical_model + File.separator +
                        Report.options_intern.f_out_validation_logos_biophysical_png);

                File tempDir = FileManagement.getFileIfInDirectory(sourceDir, "[0-9]+_" + tfGroup.getName(), false);

                File sourceFile = FileManagement.getFileIfInDirectory(tempDir, ".+\\.png$", true);

                if (sourceFile != null)
                {
                    FileManagement.copyFile(sourceFile, targetFile);

                    ArrayList<ArrayList<String>> options = new ArrayList<>();

                    options.add(new ArrayList<>(List.of(targetFile.getName())));

                    frame = frame.replace("{BIOPHYSICAL_MODEL}",
                            StructureElements.generateImageSelector("logosBiophysicalModel", targetFile.getParentFile(),
                                    options));
                }
            } // BIOPHYSICAL MODEL

            {
                File sDir = new File(sourceParentDir.getAbsolutePath() + File.separator +
                        Report.options_intern.folder_out_distribution_logos_TF_sequence);

                File targetDir = new File(Report.options_intern.com2pose_working_directory + File.separator +
                        Report.options_intern.d_out_validation + File.separator + tfGroup.getName() + File.separator +
                        Report.options_intern.d_out_validation_logos_tf_sequence);

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
                        ArrayList<ArrayList<String>> options = new ArrayList<>();
                        options.add(files);

                        frame = frame.replace("{TF_SEQUENCE}",
                                StructureElements.generateImageSelector("logosTfSequence", targetDir, options));
                    }
                }
            } // TF Sequence

            {
                File sDir = new File(sourceParentDir.getAbsolutePath() + File.separator +
                        Report.options_intern.folder_out_distribution_logos_binding_sequence);

                File targetDir = new File(Report.options_intern.com2pose_working_directory + File.separator +
                        Report.options_intern.d_out_validation + File.separator + tfGroup.getName() + File.separator +
                        Report.options_intern.d_out_validation_logos_tf_binding_sequence);

                File sourceDir = FileManagement.getFileIfInDirectory(sDir, "[0-9]+_" + tfGroup.getName(), false);

                if (sourceDir != null)
                {
                    FileManagement.copyDirectory(sourceDir, targetDir, true, ".*_biophysical_model.png$",
                            Arrays.asList(tfGroup.getName(), "biophysical_model"));

                    frame = frame.replace("{TF_BINDING_SEQUENCE}",
                            StructureElements.generateImageSelector("logosTfBindingSequence", targetDir,
                                    List.of(SelectorTypes.DISTRIBUTION_OPTIONS)));
                }
            } // TF binding sequence
        } // LOGOS

        {
            File sourceDir = new File(Report.options_intern.com2pose_working_directory + File.separator +
                    Report.options_intern.folder_out_igv + File.separator +
                    Report.options_intern.folder_out_igv_dcg_target_genes);

            File targetDir = new File(Report.options_intern.com2pose_working_directory + File.separator +
                    Report.options_intern.d_out_validation);

            File source = FileManagement.getFileIfInDirectory(sourceDir, tfGroup.getName(), false);

            String id = "validationIGV";
            File target =
                    new File(targetDir.getAbsolutePath() + File.separator + tfGroup.getName() + File.separator + id);

            frame = frame.replace("{VALIDATION_IGV_DISABLED}", (source == null) ? "disabled" : "");

            if (source != null)
            {
                FileManagement.copyDirectory(source, target, true, ".+\\.png$", new ArrayList<>());
            }
            frame = frame.replace("{VALIDATION_IGV}", (source == null) ? "" :
                    StructureElements.generateImageSelector(id, target,
                            Arrays.asList(SelectorTypes.HISTONE_MODIFICATIONS, SelectorTypes.GROUP_PAIRINGS,
                                    SelectorTypes.EMPTY_DROPDOWN)));
        } // IGV

        frame = StructureElements.setGeneCardLinks(frame, tfGroup);

        FileManagement.writeHTML(Report.options_intern.com2pose_working_directory + File.separator +
                Report.options_intern.d_out_validation + File.separator + tfGroup.getName() + File.separator +
                tfGroup.getName() + ".html", frame, 2);

        return true;
    }

    static boolean generateDistribution(TranscriptionFactorGroup tfGroup) throws IOException
    {
        String id = "distributionPlots";

        File templateFile = new File(Report.options_intern.path_to_COM2POSE + File.separator +
                Report.options_intern.f_report_resources_distribution_distribution_html);

        File d_distribution_output = new File(Report.options_intern.com2pose_working_directory + File.separator +
                Report.options_intern.d_out_distribution + File.separator + tfGroup.getName());

        File d_distribution_plots = new File(d_distribution_output.getAbsolutePath() + File.separator + id);


        File d_plots_hm = new File(Report.options_intern.com2pose_working_directory + File.separator +
                Report.options_intern.folder_out_distribution + File.separator +
                Report.options_intern.folder_out_distribution_plots + File.separator +
                Report.options_intern.folder_out_distribution_plots_HM);

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

        String frame =
                StructureElements.getFrame(tfGroup.getName() + " - Distribution", templateFile.getAbsolutePath());

        frame = frame.replace("{TFNAME}", tfGroup.getName());
        frame = StructureElements.setBasicData(frame, tfGroup);


        frame = frame.replace("{DISTRIBUTION_PLOTS}", StructureElements.generateImageSelector(id, d_distribution_plots,
                List.of(SelectorTypes.DISTRIBUTION_OPTIONS)));


        frame = StructureElements.setGeneCardLinks(frame, tfGroup);


        {
            File hmsDir = new File(Report.options_intern.com2pose_working_directory + File.separator +
                    Report.options_intern.folder_out_distribution + File.separator +
                    Report.options_intern.folder_out_distribution_stats + File.separator +
                    Report.options_intern.folder_out_distribution_stats_HM);
            HashMap<String, Integer> ranks = new HashMap<>();
            HashMap<String, Integer> sizes = new HashMap<>();

            for (File hmDir : Objects.requireNonNull(hmsDir.listFiles()))
            {
                if (!hmDir.isDirectory())
                {
                    continue;
                }

                File statsFile = new File(hmDir.getAbsolutePath() + File.separator +
                        Report.options_intern.file_suffix_distribution_analysis_plot_stats);

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
                        String[] lines = FileManagement.loadFile(statsFile.getAbsolutePath()).split("\n");
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
                score = Double.parseDouble(FileManagement.findValueInTable(tfGroup.getName(), 1, 2, new File(
                        Report.options_intern.com2pose_working_directory + File.separator +
                                Report.options_intern.folder_out_distribution + File.separator +
                                Report.options_intern.folder_out_distribution_dcg + File.separator +
                                Report.options_intern.file_suffix_distribution_analysis_dcg), "\t", true));
                frame = frame.replace("{DCG_SCORE}", Report.formatter.format(score));
            } catch (NoSuchFieldException ignored)
            {
                frame = frame.replace("{DCG_SCORE}", "-");
            }

            frame = frame.replace("{DISCOUNTED_CUMULATIVE_GAIN}", sb_dcg.toString());
        }

        FileManagement.writeHTML(d_distribution_output + File.separator + tfGroup.getName() + ".html", frame, 2);

        return true;
    }

    static boolean generateRegression(TranscriptionFactorGroup tfGroup) throws IOException
    {
        File templateFile = new File(Report.options_intern.path_to_COM2POSE + File.separator +
                Report.options_intern.f_report_resources_regression_regression_html);

        File d_in_plots = new File(
                Report.options_intern.com2pose_working_directory + File.separator + Report.options_intern.folder_plots);

        File d_out_regression = new File(Report.options_intern.com2pose_working_directory + File.separator +
                Report.options_intern.d_out_regression);

        String frame = StructureElements.getFrame(tfGroup.getName() + " - Regression", templateFile.getAbsolutePath());

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

        FileManagement.writeHTML(Report.options_intern.com2pose_working_directory + File.separator +
                Report.options_intern.d_out_regression + File.separator + tfGroup.getName() + File.separator +
                tfGroup.getName() + ".html", frame, 2);

        return true;
    }
}
