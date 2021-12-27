package util;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.*;

import static java.nio.file.StandardCopyOption.REPLACE_EXISTING;

public class Report
{
    private final Logger logger;
    private final Options_intern options_intern;
    private final ArrayList<TranscriptionFactorGroup> transcriptionFactorGroups = new ArrayList<>();
    private final DecimalFormat formatter = new DecimalFormat("0.###");

    public Report(Options_intern options_intern) throws IOException
    {
        this.options_intern = options_intern;
        logger = new Logger(true, options_intern.com2pose_working_directory);
        loadTFs();
    }

    private record TranscriptionFactor(String geneID, String name, Map<String, Map<String, Number>> log2fc,
                                       Map<String, Number> tpm, Map<String, Number> normex,
                                       ArrayList<String> histoneModifications, boolean hasIGV)
    {
    }

    private record TranscriptionFactorGroup(String name, ArrayList<TranscriptionFactor> transcriptionFactors,
                                            boolean realGroup, boolean hasIGV)
    {
    }

    private void loadTFs() throws FileNotFoundException
    {
        File tf_file = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_distribution +
                        File.separator + options_intern.folder_out_distribution_dcg + File.separator +
                        options_intern.file_suffix_distribution_analysis_dcg);

        File geneIDFile = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_deseq2_preprocessing + File.separator +
                options_intern.file_suffix_deseq2_mapping);

        File d_plots =
                new File(options_intern.com2pose_working_directory + File.separator + options_intern.folder_plots);

        ArrayList<String> histoneModifications = new ArrayList<>();

        for (File f_histoneModification : Objects.requireNonNull(d_plots.listFiles()))
        {
            histoneModifications.add(f_histoneModification.getName());
        }

        try (Scanner scanner = new Scanner(tf_file))
        {
            boolean firstLine = true;
            while (scanner.hasNextLine())
            {
                ArrayList<TranscriptionFactor> tf_group = new ArrayList<>();
                String line = scanner.nextLine();
                if (firstLine)
                {
                    firstLine = false;
                    continue;
                }
                String tfGroupName = line.split("\t")[1];

                for (String tf_name : tfGroupName.split("\\.\\."))
                {
                    try
                    {
                        String geneID = findValueInTable(tf_name, 1, 0, geneIDFile, "\t", true);
                        Map<String, Map<String, Number>> log2fc = new HashMap<>();
                        Map<String, Number> tpm = new HashMap<>();
                        Map<String, Number> normex = new HashMap<>();
                        boolean hasIGV = false;

                        {   //LOG2FC
                            File d_log2fs = new File(options_intern.com2pose_working_directory + File.separator +
                                    options_intern.folder_name_deseq2_output);

                            for (File entry : Objects.requireNonNull(d_log2fs.listFiles()))
                            {
                                if (entry.isFile())
                                {
                                    String group1 = entry.getName().split("_")[0];
                                    String group2 = entry.getName().split("_")[1];

                                    try
                                    {
                                        double log2fc_value =
                                                Double.parseDouble(findValueInTable(geneID, 0, 1, entry, "\t", false));

                                        if (!log2fc.containsKey(group1))
                                        {
                                            log2fc.put(group1, new HashMap<>());
                                        }
                                        if (!log2fc.containsKey(group2))
                                        {
                                            log2fc.put(group2, new HashMap<>());
                                        }

                                        log2fc.get(group1).put(group2, log2fc_value);
                                        log2fc.get(group2).put(group1, log2fc_value);
                                    } catch (NoSuchFieldException ignored)
                                    {
                                    }
                                }
                            }
                        }   //LOG2FC

                        {   //TPM
                            File d_tpm = new File(options_intern.com2pose_working_directory + File.separator +
                                    options_intern.folder_name_deseq2_preprocessing + File.separator +
                                    options_intern.folder_name_deseq2_preprocessing_tpm + File.separator +
                                    options_intern.folder_name_deseq2_preprocessing_tpm_results);

                            for (File entry : Objects.requireNonNull(d_tpm.listFiles()))
                            {
                                if (entry.isFile() && entry.getName().endsWith(".csv"))
                                {
                                    String group = entry.getName().split("_")[0];

                                    try
                                    {
                                        double tpm_value =
                                                Double.parseDouble(findValueInTable(geneID, 0, 3, entry, "\t", false));
                                        tpm.put(group, tpm_value);
                                    } catch (NoSuchFieldException ignored)
                                    {
                                    }
                                }
                            }
                        }   //TPM

                        {   //Normalized expression
                            File d_normex = new File(options_intern.com2pose_working_directory + File.separator +
                                    options_intern.folder_name_deseq2_preprocessing + File.separator +
                                    options_intern.folder_name_deseq2_preprocessing_gene_symbols);


                            for (File entry : Objects.requireNonNull(d_normex.listFiles()))
                            {
                                if (entry.isFile() && entry.getName().endsWith(".csv"))
                                {
                                    String group = entry.getName().split("\\.")[0];

                                    try
                                    {
                                        int exp_value =
                                                Integer.parseInt(findValueInTable(geneID, 1, 2, entry, "\t", false));
                                        normex.put(group, exp_value);
                                    } catch (NoSuchFieldException ignored)
                                    {
                                    }
                                }
                            }
                        }   //Normalized expression

                        {   //Check if TF has IGV
                            File d_igv_screenshots = new File(
                                    options_intern.com2pose_working_directory + File.separator +
                                            options_intern.folder_out_igv);

                            if (d_igv_screenshots.isDirectory())
                            {
                                File d_tf_igv = new File(
                                        d_igv_screenshots + File.separator + options_intern.folder_out_igv_own_data);

                                for (File entry : Objects.requireNonNull(d_tf_igv.listFiles()))
                                {
                                    String name = entry.getName().split("_")[1];
                                    if (name.equals(tf_name))
                                    {
                                        hasIGV = true;
                                        break;
                                    }
                                }
                            }
                        }   //Check if TF has IGV

                        TranscriptionFactor tf =
                                new TranscriptionFactor(geneID, tf_name, log2fc, tpm, normex, histoneModifications,
                                        hasIGV);
                        tf_group.add(tf);
                    } catch (NoSuchFieldException ignore)
                    {
                        System.out.println("Not found: " + tf_name);
                    }
                }

                boolean groupHasIGV = false;

                if (tf_group.size() > 1)
                {
                    {   //Check if TF has IGV
                        File d_igv_screenshots = new File(options_intern.com2pose_working_directory + File.separator +
                                options_intern.folder_out_igv);

                        if (d_igv_screenshots.isDirectory())
                        {
                            File d_tf_igv = new File(
                                    d_igv_screenshots + File.separator + options_intern.folder_out_igv_own_data);

                            for (File entry : Objects.requireNonNull(d_tf_igv.listFiles()))
                            {
                                String name = entry.getName().split("_")[1];
                                if (name.equals(tfGroupName))
                                {
                                    groupHasIGV = true;
                                    break;
                                }
                            }
                        }
                    }   //Check if TF has IGV
                }

                if (tf_group.size() > 0)
                {
                    transcriptionFactorGroups.add(
                            new TranscriptionFactorGroup(tfGroupName, tf_group, tf_group.size() > 1, groupHasIGV));
                }
            }
        }
    }

    public void generate() throws IOException
    {
        logger.logLine("[REPORT] Start generating report");

        generateHome();
        styleAndScript();
        generateParameters();

        for (TranscriptionFactorGroup tfGroup : transcriptionFactorGroups)
        {
            generateValidation(tfGroup);
            generateDistribution(tfGroup);
            generateRegression(tfGroup);
        }

        logger.logLine("[REPORT] Finished generating report");
    }

    private void generateParameters() throws IOException
    {
        logger.logLine("[REPORT] Start generating report parameters page");
        String parameters = loadFrame();
        parameters = parameters.replace("{BODY}", loadFile(options_intern.path_to_COM2POSE + File.separator +
                options_intern.f_report_resources_parameters_parameters_html));
        parameters = parameters.replace("{TITLE}", "Parameters");
        parameters = parameters.replace("{TOOLS}", loadFile(options_intern.path_to_COM2POSE + File.separator +
                options_intern.f_report_resources_parameters_tool_html));
        parameters = parameters.replace("{PARAMETERS}", loadFile(options_intern.path_to_COM2POSE + File.separator +
                options_intern.f_report_resources_parameters_parameter_html));

        parameters = relativate(parameters, 0);

        writeFile(options_intern.com2pose_working_directory + File.separator + options_intern.f_out_report_parameters,
                parameters);
        logger.logLine("[REPORT] Finished generating report parameters page");
    }

    private void generateHome() throws IOException
    {
        logger.logLine("[REPORT] Start generating report home page");
        String home = loadFrame();
        home = home.replace("{BODY}", loadFile(
                options_intern.path_to_COM2POSE + File.separator + options_intern.f_report_resources_home_home_html));
        home = home.replace("{TITLE}", "Overview");


        StringBuilder sb_tfs = new StringBuilder();


        int i = 1;

        for (TranscriptionFactorGroup tfGroup : transcriptionFactorGroups)
        {
            ArrayList<TranscriptionFactor> tfList = tfGroup.transcriptionFactors;

            if (!tfGroup.realGroup)
            {
                for (TranscriptionFactor transcriptionFactor : tfList)
                {
                    String tf_string = loadFile(options_intern.path_to_COM2POSE + File.separator +
                            options_intern.f_report_resources_home_tf_html);

                    tf_string = tf_string.replace("{BUTTONBAR}", getButtonBar(transcriptionFactor));

                    tf_string = tf_string.replace("{TF_NAME}", i + ". " + transcriptionFactor.name);

                    tf_string = tf_string.replace("{BASICDATA}", getBasicData(transcriptionFactor));

                    tf_string = tf_string.replace("{GENEID}", transcriptionFactor.geneID);

                    tf_string = tf_string.replace("{HASVALIDATION}", transcriptionFactor.hasIGV ? "" : "disabled");
                    tf_string = tf_string.replace("{HASDISTRIBUTION}", "disabled");
                    tf_string = tf_string.replace("{HASREGRESSION}", "");

                    sb_tfs.append(tf_string);
                }
            } else
            {
                String tfGroupString = loadFile(options_intern.path_to_COM2POSE + File.separator +
                        options_intern.f_report_resources_home_tfGroup_html);

                tfGroupString = tfGroupString.replace("{TF_NAME}", i + ". " + tfGroup.name);
                tfGroupString = tfGroupString.replace("{BUTTONBAR}", getButtonBar(tfGroup));
                tfGroupString = tfGroupString.replace("{SINGLE_TFS}", getBasicData(tfGroup));

                sb_tfs.append(tfGroupString);
            }

            i++;
        }
        home = home.replace("{TFS}", sb_tfs.toString());

        home = relativate(home, 0);


        writeFile(options_intern.com2pose_working_directory + File.separator + options_intern.f_out_report_home, home);
        logger.logLine("[REPORT] Finished generating report home page");
    }

    private String getBasicDataEntry(String name, Map<String, Number> data) throws IOException
    {
        if (data.size() == 0)
        {
            return "";
        }

        String template = loadFile(options_intern.path_to_COM2POSE + File.separator +
                options_intern.f_report_resources_basicdata_entry_html);

        template = template.replace("{NAME}", name);

        StringBuilder sb_data = new StringBuilder();

        sb_data.append("<div class='keyvaluepaircontainer'>");

        for (Map.Entry<String, Number> kvPair : data.entrySet())
        {
            sb_data.append("<div class=\"keyvaluepair\"><h4>" + kvPair.getKey() + "</h4><p>" +
                    formatter.format(kvPair.getValue()) + "</p></div>");
        }

        sb_data.append("</div>");

        template = template.replace("{DATA}", sb_data.toString());

        return template;
    }

    private String getBasicDataLog2fc(String tfName, Map<String, Map<String, Number>> log2fc) throws IOException
    {
        if (log2fc.size() == 0)
        {
            return "";
        }

        String template = loadFile(options_intern.path_to_COM2POSE + File.separator +
                options_intern.f_report_resources_basicdata_entry_html);

        template = template.replace("{NAME}", "LOG2FC");

        StringBuilder sb_data = new StringBuilder();

        List<String> groups = new ArrayList<>(log2fc.keySet());
        Collections.sort(groups);

        sb_data.append("<table>");
        sb_data.append("<tr>");
        sb_data.append("<th></th>");
        int i = 0;
        for (String group : groups)
        {
            sb_data.append("<th id='").append(tfName).append("-col-").append(i).append("'>");
            sb_data.append(group);
            sb_data.append("</th>");
            i++;
        }
        sb_data.append("</tr>");

        i = 0;
        for (String group : groups)
        {
            sb_data.append("<tr>");

            sb_data.append("<th id='").append(tfName).append("-row-").append(i).append("'>");
            sb_data.append(group);
            sb_data.append("</th>");

            int j = 0;
            for (String match : groups)
            {
                String parameters = "\"" + tfName + "\", " + j + ", " + i;
                sb_data.append("<td onmouseover='tableMouseOver(").append(parameters)
                        .append(")' onmouseout" + "='tableMouseOut(").append(parameters).append(")'>");
                if (!match.equals(group))
                {
                    sb_data.append(formatter.format(log2fc.get(group).get(match)));
                }
                sb_data.append("</td>");

                j++;
            }

            sb_data.append("</tr>");

            i++;
        }

        sb_data.append("</table>");

        template = template.replace("{DATA}", sb_data.toString());

        return template;
    }

    private String getBasicData(TranscriptionFactor transcriptionFactor) throws IOException
    {
        String template = loadFile(
                options_intern.path_to_COM2POSE + File.separator + options_intern.f_report_resources_basicdata_html);


        template =
                template.replace("{LOG2FC}", getBasicDataLog2fc(transcriptionFactor.name, transcriptionFactor.log2fc));

        template = template.replace("{TPM}", getBasicDataEntry("TPM", transcriptionFactor.tpm));

        template = template.replace("{NORMEX}", getBasicDataEntry("Norm. expression", transcriptionFactor.normex));

        return template;
    }

    private String getBasicData(TranscriptionFactorGroup tfGroup) throws IOException
    {
        StringBuilder basicData = new StringBuilder();
        String tfTemplate = loadFile(
                options_intern.path_to_COM2POSE + File.separator + options_intern.f_report_resources_home_tf_html);

        for (TranscriptionFactor tf : tfGroup.transcriptionFactors)
        {
            String tfString = tfTemplate.replace("{BASICDATA}", getBasicData(tf));
            tfString = tfString.replace("{GENEID}", tf.geneID);
            tfString = tfString.replace("{BUTTONBAR}", "");
            tfString = tfString.replace("{TF_NAME}", tf.name);
            basicData.append(tfString);
        }

        return basicData.toString();
    }

    private String getButtonBar(TranscriptionFactor tf) throws IOException
    {
        return getButtonBar(tf.name, tf.hasIGV, true, true);
    }

    private String getButtonBar(TranscriptionFactorGroup tfGroup) throws IOException
    {
        return getButtonBar(tfGroup.name, tfGroup.hasIGV, true, true);
    }

    private String getButtonBar(String name, boolean hasValidation, boolean hasDistribution, boolean hasRegression)
            throws IOException
    {
        String buttonbar = loadFile(options_intern.path_to_COM2POSE + File.separator +
                options_intern.f_report_resources_home_buttonbar_html);

        buttonbar = buttonbar.replace("{VALIDATION}",
                "VALIDATION" + File.separator + name + File.separator + name + ".html");
        buttonbar = buttonbar.replace("{DISTRIBUTION}",
                "DISTRIBUTION" + File.separator + name + File.separator + name + ".html");
        buttonbar = buttonbar.replace("{REGRESSION}",
                "REGRESSION" + File.separator + name + File.separator + name + ".html");

        buttonbar = buttonbar.replace("{HASVALIDATION}", hasValidation ? "" : "disabled");
        buttonbar = buttonbar.replace("{HASREGRESSION}", hasRegression ? "" : "disabled");
        buttonbar = buttonbar.replace("{HASDISTRIBUTION}", hasDistribution ? "" : "disabled");

        return buttonbar;
    }


    private void generateValidation(TranscriptionFactorGroup tfGroup) throws IOException
    {
        if (tfGroup.realGroup)
        {
            generateValidation(tfGroup.name, getBasicData(tfGroup));
        } else
        {
            TranscriptionFactor tf = tfGroup.transcriptionFactors.get(0);

            generateValidation(tf.name, getBasicData(tf));
        }
    }

    private void generateValidation(String name, String basicData) throws IOException
    {
        File templateFile = new File(options_intern.path_to_COM2POSE + File.separator +
                options_intern.f_report_resources_validation_validation_html);

        File d_igv_screenshots = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_igv +
                        File.separator + options_intern.folder_out_igv_own_data);

        File d_out_validation =
                new File(options_intern.com2pose_working_directory + File.separator + options_intern.d_out_validation);

        String frame = loadFrame();
        frame = frame.replace("{BODY}", loadFile(templateFile.getAbsolutePath()));

        frame = frame.replace("{TFNAME}", name);

        frame = frame.replace("{TITLE}", name + " - Validation");

        frame = frame.replace("{BASICDATA}", basicData);

        File source = null;

        for (File entry : d_igv_screenshots.listFiles())
        {
            String entryName = entry.getName();
            entryName = entryName.split("_")[1];
            if (entryName.equals(name))
            {
                source = entry;
                break;
            }
        }

        File target = new File(d_out_validation.getAbsolutePath() + File.separator + name);

        if (!(source == null))
        {
            String three_level_image_selector =
                    generateThreeLevelImageSelector("validationPlot", source, target, new ArrayList<>(List.of(name)));

            frame = frame.replace("{VALIDATION_OWN_TF}", three_level_image_selector);

            frame = relativate(frame, 2);

            frame = frame.replace("{GENECARD}", options_intern.link_report_genecards.replace("{GENE}", name));

            writeFile(options_intern.com2pose_working_directory + File.separator + options_intern.d_out_validation +
                    File.separator + name + File.separator + name + ".html", frame);
        }
    }

    private String generateThreeLevelImageSelector(String name, File sourceDir, File targetDir) throws IOException
    {
        return generateThreeLevelImageSelector(name, sourceDir, targetDir, new ArrayList<>());
    }

    private String generateThreeLevelImageSelector(String name, File sourceDir, File targetDir,
                                                   ArrayList<String> specialRemovables) throws IOException
    {
        String suffix = ".png";
        HashMap<String, HashMap<String, ArrayList<String>>> combinations = new HashMap<>();

        for (File group_dir : Objects.requireNonNull(sourceDir.listFiles()))
        {
            if (group_dir.isFile())
            {
                continue;
            }

            String group = group_dir.getName();

            if (group.equals("A_SESSIONS"))
            {
                continue;
            }

            combinations.put(group, new HashMap<>());

            for (File subgroup_dir : Objects.requireNonNull(group_dir.listFiles()))
            {
                if (subgroup_dir.isFile())
                {
                    continue;
                }

                String subgroup = subgroup_dir.getName();

                if (subgroup.split("_").length > 1)
                {
                    subgroup = subgroup.split("_")[1];
                }

                combinations.get(group).put(subgroup, new ArrayList<>());

                for (File image_file : Objects.requireNonNull(subgroup_dir.listFiles()))
                {
                    ArrayList<String> removables = new ArrayList<>(List.of(group, subgroup, suffix, "threshold"));
                    removables.addAll(specialRemovables);

                    if (image_file.isDirectory() || !image_file.getName().endsWith(suffix))
                    {
                        continue;
                    }

                    String relevantFileName = image_file.getName();

                    for (String entry : removables)
                    {
                        relevantFileName = relevantFileName.replace(entry, "");
                    }

                    if (relevantFileName.matches("[0-9]+_.*_.*")) // Check for validation only
                    {
                        relevantFileName = relevantFileName.replaceAll("[0-9]+_", "");
                    }

                    while (relevantFileName.endsWith("_"))
                    {
                        relevantFileName = relevantFileName.substring(0, relevantFileName.length() - 1);
                    }
                    while (relevantFileName.startsWith("_"))
                    {
                        relevantFileName = relevantFileName.substring(1);
                    }

                    copyFile(image_file, new File(
                            targetDir.getAbsolutePath() + File.separator + group + File.separator + subgroup +
                                    File.separator + relevantFileName + suffix));

                    combinations.get(group).get(subgroup).add(relevantFileName);
                }
            }
        }

        String three_level_image_selector = loadFile(options_intern.f_report_resources_three_level_image_selector_html);

        Set<String> groups = combinations.keySet();
        Set<String> subgroups = new HashSet<>();

        StringBuilder sb_groups = new StringBuilder();

        for (String group : groups)
        {
            subgroups.addAll(combinations.get(group).keySet());

            sb_groups.append("<button class=\"{ID} group-selector\" onclick=\"select_group" + "('{ID}', this, " +
                    "{ID}Combinations)\" " + "value=\"" + group + "\">" + group + "</button>");
        }
        three_level_image_selector = three_level_image_selector.replace("{GROUPS}", sb_groups.toString());

        StringBuilder sb_subgroups = new StringBuilder();
        for (String subgroup : subgroups)
        {
            sb_subgroups.append(
                    "<button class=\"{ID} subgroup-selector\" onclick=\"select_subgroup" + "('{ID}', " + "this, " +
                            "{ID}Combinations)\" " + "value=\"" + subgroup + "\">" + subgroup + "</button>");
        }

        three_level_image_selector = three_level_image_selector.replace("{SUBGROUPS}", sb_subgroups.toString());

        String json;
        {
            HashMap<String, HashMap<String, String>> lv1 = new HashMap<>();
            HashMap<String, String> lv2 = new HashMap<>();

            for (String hm : combinations.keySet())
            {
                lv1.put(hm, new HashMap<>());
                for (String group : combinations.get(hm).keySet())
                {
                    StringBuilder sb_genes = new StringBuilder("[");
                    for (String gene : combinations.get(hm).get(group))
                    {
                        sb_genes.append("\"");
                        sb_genes.append(gene);
                        sb_genes.append("\",");
                    }
                    sb_genes.setLength(sb_genes.length() - 1);
                    sb_genes.append("]");
                    lv1.get(hm).put(group, sb_genes.toString());
                }
                lv2.put(hm, mapToJson(lv1.get(hm)));
            }

            json = mapToJson(lv2);
        }

        three_level_image_selector = three_level_image_selector.replace("{ID}", name);

        three_level_image_selector = three_level_image_selector.replace("{COMBINATIONS}", json);

        return three_level_image_selector;
    }

    private void generateDistribution(TranscriptionFactorGroup tfGroup) throws IOException
    {
        if (tfGroup.realGroup)
        {
            generateDistribution(tfGroup.name, getBasicData(tfGroup));
        } else
        {
            TranscriptionFactor tf = tfGroup.transcriptionFactors.get(0);

            generateDistribution(tf.name, getBasicData(tf));
        }
    }

    private void generateDistribution(String name, String basicData) throws IOException
    {
        File templateFile = new File(options_intern.path_to_COM2POSE + File.separator +
                options_intern.f_report_resources_distribution_distribution_html);

        File d_distribution_output = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.d_out_distribution +
                        File.separator + name);

        ArrayList<String> existingHMs = new ArrayList<>();

        File d_plots = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_distribution +
                        File.separator + options_intern.folder_out_distribution_plots + File.separator +
                        options_intern.folder_out_distribution_plots_HM);

        for (File d_hm : d_plots.listFiles())
        {
            for (File f_plot : d_hm.listFiles())
            {
                if (f_plot.getName().split("\\.")[0].equals(name))
                {
                    existingHMs.add(d_hm.getName());
                    copyFile(f_plot, new File(
                            d_distribution_output.getAbsolutePath() + File.separator + d_hm.getName() + ".png"));
                    break;
                }
            }
        }

        if (existingHMs.size() == 0)
        {
            return;
        }

        String frame = loadFrame();
        frame = frame.replace("{BODY}", loadFile(templateFile.getAbsolutePath()));

        frame = frame.replace("{TFNAME}", name);
        frame = frame.replace("{BASICDATA}", basicData);

        frame = frame.replace("{TITLE}", name + " - Distribution");

        StringBuilder sb_histoneModifications = new StringBuilder();

        boolean first = true;

        for (String histoneModification : existingHMs)
        {
            sb_histoneModifications.append(
                    "<button class=\"selector" + (first ? " active" : "") + "\" name=\"distribution-plot\" value=\"" +
                            histoneModification + "\">" + histoneModification + "</button>");
            first = false;
        }

        frame = frame.replace("{HISTONEMODIFICATIONS}", sb_histoneModifications.toString());
        frame = frame.replace("{DISTRIBUTION-PLOT-PATH}", existingHMs.get(0) + ".png");

        frame = frame.replace("{GENECARD}", options_intern.link_report_genecards.replace("{GENE}", name));

        frame = relativate(frame, 2);

        writeFile(d_distribution_output + File.separator + name + ".html", frame);
    }

    private void generateRegression(TranscriptionFactorGroup tfGroup) throws IOException
    {
        if (tfGroup.realGroup)
        {
            generateRegression(tfGroup.name);
        } else
        {
            TranscriptionFactor tf = tfGroup.transcriptionFactors.get(0);

            generateRegression(tf.name);
        }
    }

    private void generateRegression(String name) throws IOException
    {

        File templateFile = new File(options_intern.path_to_COM2POSE + File.separator +
                options_intern.f_report_resources_regression_regression_html);

        File d_in_plots =
                new File(options_intern.com2pose_working_directory + File.separator + options_intern.folder_plots);

        File d_out_regression =
                new File(options_intern.com2pose_working_directory + File.separator + options_intern.d_out_regression);

        String frame = loadFrame();

        frame = frame.replace("{BODY}", loadFile(templateFile.getAbsolutePath()));

        frame = frame.replace("{TFNAME}", name);

        frame = frame.replace("{TITLE}", name + " - Regression");

        File target = new File(d_out_regression.getAbsolutePath() + File.separator + name);
        String three_level_image_selector = generateThreeLevelImageSelector("regressionPlot", d_in_plots, target);

        frame = frame.replace("{HEATMAPS}", three_level_image_selector);

        frame = frame.replace("{ADDHEAD}", "<script src=\"COMBINATIONS.js\"></script>");

        frame = relativate(frame, 2);

        writeFile(options_intern.com2pose_working_directory + File.separator + options_intern.d_out_regression +
                File.separator + name + File.separator + name + ".html", frame);
    }

    private String findValueInTable(String term, int searchIndex, int resultIndex, File file, String sep,
                                    boolean ignoreCase) throws FileNotFoundException, NoSuchFieldException
    {
        if (file.isFile())
        {
            try (Scanner scanner = new Scanner(file))
            {
                while (scanner.hasNextLine())
                {
                    String line = scanner.nextLine();
                    if (line.split(sep).length > searchIndex && line.split(sep).length > resultIndex)
                    {

                        if (term.equals(line.split(sep)[searchIndex]) ||
                                ignoreCase && term.equalsIgnoreCase(line.split(sep)[searchIndex]))
                        {
                            return line.split(sep)[resultIndex];
                        }
                    }
                }
            }
        }
        throw new NoSuchFieldException();
    }

    private String relativate(String content, int depth)
    {
        content = content.replace("{RELATIVATION}", (".." + File.separator).repeat(depth));
        return content;
    }

    private void styleAndScript() throws IOException
    {
        logger.logLine("[REPORT] Start moving CSS and JS components");
        String css =
                loadFile(options_intern.path_to_COM2POSE + File.separator + options_intern.f_report_resources_style);
        String script =
                loadFile(options_intern.path_to_COM2POSE + File.separator + options_intern.f_report_resources_script);

        writeFile(options_intern.com2pose_working_directory + File.separator + options_intern.f_out_report_style, css);
        writeFile(options_intern.com2pose_working_directory + File.separator + options_intern.f_out_report_script,
                script);
    }

    private void copyFile(File source, File target) throws IOException
    {
        if (target.exists())
        {
            return;
        }

        if (source.exists())
        {
            if (!target.getParentFile().exists())
            {
                target.getParentFile().mkdirs();
            }

            Files.copy(source.toPath(), target.toPath(), REPLACE_EXISTING);
        }
    }

    private void writeFile(String path, String content) throws IOException
    {
        content = content.replace("{ADDHEAD}", "");
        File target = new File(path);
        if (!target.exists())
        {
            target.getParentFile().mkdirs();
            target.createNewFile();
        }
        Files.writeString(target.toPath(), content);
    }

    private String loadFrame() throws IOException
    {
        return loadFile(
                options_intern.path_to_COM2POSE + File.separator + options_intern.f_report_resources_frame_html);
    }

    private String loadFile(String path) throws IOException
    {
        File source = new File(path);
        return Files.readString(source.toPath());
    }

    private String mapToJson(Map<String, String> map)
    {
        StringBuilder sb_output = new StringBuilder("{");

        for (Map.Entry<String, String> entry : map.entrySet())
        {
            boolean valueIsJson = (entry.getValue().startsWith("{") && entry.getValue().endsWith("}")) ||
                    (entry.getValue().startsWith("[") && entry.getValue().endsWith("]"));
            sb_output.append("\"");
            sb_output.append(entry.getKey());
            sb_output.append("\":");
            if (!valueIsJson)
            {
                sb_output.append("\"");
            }
            sb_output.append(entry.getValue());
            if (!valueIsJson)
            {
                sb_output.append("\"");
            }
            sb_output.append(",");
        }
        sb_output.setLength(sb_output.length() - 1);
        sb_output.append("}");

        return sb_output.toString();
    }
}
