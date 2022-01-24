package util;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.*;
import java.util.concurrent.*;

import static java.nio.file.StandardCopyOption.REPLACE_EXISTING;

public class Report
{
    private final Logger logger;
    private final Options_intern options_intern;
    private final ArrayList<TranscriptionFactorGroup> transcriptionFactorGroups = new ArrayList<>();
    private final DecimalFormat formatter = new DecimalFormat("0.###");
    private final ExecutorService executorService = Executors.newFixedThreadPool(10);

    public Report(Options_intern options_intern) throws IOException
    {
        this.options_intern = options_intern;
        logger = new Logger(true, options_intern.com2pose_working_directory);
        loadTFs();
    }

    private record TranscriptionFactor(String geneID, String name, Map<String, Map<String, Number>> log2fc,
                                       Map<String, Number> tpm, Map<String, Number> normex,
                                       ArrayList<String> histoneModifications)
    {
    }

    private static class TranscriptionFactorGroup
    {
        boolean hasValidation, hasDistribution, hasRegression;
        private final String name;
        private final ArrayList<TranscriptionFactor> transcriptionFactors;
        Map<String, Map<String, Number>> regressionCoefficients;
        boolean realGroup;

        TranscriptionFactorGroup(String name, ArrayList<TranscriptionFactor> transcriptionFactors,
                                 Map<String, Map<String, Number>> regressionCoefficients, boolean realGroup)
        {
            this.name = name;
            this.transcriptionFactors = transcriptionFactors;
            this.regressionCoefficients = regressionCoefficients;
            this.realGroup = realGroup;
        }
    }

    private Map<String, Map<String, Map<String, Double>>> loadRegressionCoefficients() throws IOException
    {
        Map<String, Map<String, Map<String, Double>>> coefficients = new HashMap<>();

        File parentDir = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_put_DYNAMITE);

        for (File hm_dir : Objects.requireNonNull(parentDir.listFiles()))
        {
            coefficients.put(hm_dir.getName(), new HashMap<>());

            for (File groups_dir : Objects.requireNonNull(hm_dir.listFiles()))
            {
                coefficients.get(hm_dir.getName()).put(groups_dir.getName(), new HashMap<>());

                File f_data = new File(
                        groups_dir + File.separator + options_intern.file_suffix_dynamite_output_to_be_plotted);

                String data = loadFile(f_data.getAbsolutePath());

                boolean first = true;
                for (String line : data.split("\n"))
                {
                    if (first)
                    {
                        first = false;
                        continue;
                    }
                    String tf = line.split("\t")[0];
                    double value = Double.parseDouble(line.split("\t")[1]);

                    coefficients.get(hm_dir.getName()).get(groups_dir.getName()).put(tf.toUpperCase(), value);
                }
            }
        }

        return coefficients;
    }

    private void loadTFs() throws IOException
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

        Map<String, Map<String, Map<String, Double>>> allRegressionCoefficients = loadRegressionCoefficients();

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

                        TranscriptionFactor tf =
                                new TranscriptionFactor(geneID, tf_name, log2fc, tpm, normex, histoneModifications);
                        tf_group.add(tf);
                    } catch (NoSuchFieldException ignore)
                    {
                        System.out.println("Not found: " + tf_name);
                    }
                }

                Map<String, Map<String, Number>> regressionCoefficients = new HashMap<>();

                for (String hm : allRegressionCoefficients.keySet())
                {
                    regressionCoefficients.put(hm, new HashMap<>());

                    for (String group : allRegressionCoefficients.get(hm).keySet())
                    {
                        if (allRegressionCoefficients.get(hm).get(group).containsKey(tfGroupName.toUpperCase()))
                        {
                            regressionCoefficients.get(hm).put(group,
                                    allRegressionCoefficients.get(hm).get(group).get(tfGroupName.toUpperCase()));
                        }
                    }
                }

                if (tf_group.size() > 0)
                {
                    transcriptionFactorGroups.add(
                            new TranscriptionFactorGroup(tfGroupName, tf_group, regressionCoefficients,
                                    tf_group.size() > 1));
                }
            }
        }
    }

    public void generate() throws IOException, InterruptedException
    {
        logger.logLine("[REPORT] Start generating report");

        styleAndScript();
        generateParameters();

        int i = 1;

        for (TranscriptionFactorGroup tfGroup : transcriptionFactorGroups)
        {
            System.out.print("Generating report for tf " + i + " out of " + transcriptionFactorGroups.size() + ": " +
                    tfGroup.name + ".\r");
            tfGroup.hasValidation = generateValidation(tfGroup);
            tfGroup.hasDistribution = generateDistribution(tfGroup);
            tfGroup.hasRegression = generateRegression(tfGroup);
            i++;
        }

        generateHome();
        generateImportantLoci();
        generateTopLog2fc();
        generateCoOccurrence();

        executorService.shutdown();
        ThreadPoolExecutor tpe = (ThreadPoolExecutor) executorService;
        long startTime = System.currentTimeMillis();
        long totalTasks = tpe.getQueue().stream().filter(t -> !((FutureTask<?>) t).isDone()).count();

        while (!executorService.isTerminated())
        {
            long now = System.currentTimeMillis();
            long pendingTasks = tpe.getQueue().stream().filter(t -> !((FutureTask<?>) t).isDone()).count();
            double timeDelta = (now - startTime) / 1000.;
            long finishedTasks = totalTasks - pendingTasks;

            double tasksPerSecond = finishedTasks / timeDelta;
            long secondsLeft = (long) (pendingTasks / tasksPerSecond);
            long minutesLeft = secondsLeft / 60;
            secondsLeft = secondsLeft % 60;

            System.out.print(
                    "Files left to copy: " + pendingTasks + "\tETA: " + minutesLeft + "m " + secondsLeft + "s" + "\r");


            Thread.sleep(500);
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
            if (!tfGroup.realGroup)
            {
                for (TranscriptionFactor transcriptionFactor : tfGroup.transcriptionFactors)
                {
                    String tf_string = loadFile(options_intern.path_to_COM2POSE + File.separator +
                            options_intern.f_report_resources_home_tf_html);

                    tf_string = tf_string.replace("{BUTTONBAR}", getButtonBar(tfGroup));

                    tf_string = tf_string.replace("{TF_NAME}", i + ". " + transcriptionFactor.name);

                    tf_string = tf_string.replace("{BASICDATA}", getBasicData(transcriptionFactor));

                    tf_string = tf_string.replace("{GENEID}", transcriptionFactor.geneID);

                    tf_string = tf_string.replace("{HASVALIDATION}", tfGroup.hasValidation ? "" : "disabled");
                    tf_string = tf_string.replace("{HASDISTRIBUTION}", tfGroup.hasDistribution ? "" : "disabled");
                    tf_string = tf_string.replace("{HASREGRESSION}", tfGroup.hasRegression ? "" : "disabled");

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
            sb_data.append("<div class=\"keyvaluepair\"><h4>").append(kvPair.getKey()).append("</h4><p>")
                    .append(formatter.format(kvPair.getValue())).append("</p></div>");
        }

        sb_data.append("</div>");

        template = template.replace("{DATA}", sb_data.toString());

        return template;
    }

    private String getTabularData(String id, Map<String, Map<String, Number>> data)
    {
        if (data.size() == 0)
        {
            return "";
        }

        StringBuilder sb_data = new StringBuilder();

        List<String> columns = new ArrayList<>(data.keySet());
        Set<String> rowsSet = new HashSet<>();
        Collections.sort(columns);

        sb_data.append("<table>");
        sb_data.append("<tr>");
        sb_data.append("<th></th>");
        int i = 0;
        for (String column : columns)
        {
            sb_data.append("<th id='").append(id).append("-col-").append(i).append("'>");
            sb_data.append(column);
            sb_data.append("</th>");
            i++;

            rowsSet.addAll(data.get(column).keySet());
        }
        sb_data.append("</tr>");

        List<String> rows = new ArrayList<>(rowsSet);
        Collections.sort(rows);

        i = 0;
        for (String row : rows)
        {
            sb_data.append("<tr>");

            sb_data.append("<th id='").append(id).append("-row-").append(i).append("'>");
            sb_data.append(row);
            sb_data.append("</th>");

            int j = 0;
            for (String column : columns)
            {
                String parameters = "\"" + id + "\", " + j + ", " + i;
                sb_data.append("<td onmouseover='tableMouseOver(").append(parameters)
                        .append(")' onmouseout" + "='tableMouseOut(").append(parameters).append(")'>");
                if (!column.equals(row) && data.get(column).get(row) != null)
                {
                    sb_data.append(formatter.format(data.get(column).get(row)));
                } else
                {
                    sb_data.append("-");
                }
                sb_data.append("</td>");

                j++;
            }

            sb_data.append("</tr>");

            i++;
        }

        sb_data.append("</table>");

        return sb_data.toString();
    }

    private String getBasicData(TranscriptionFactor transcriptionFactor) throws IOException
    {
        String template = loadFile(
                options_intern.path_to_COM2POSE + File.separator + options_intern.f_report_resources_basicdata_html);


        String log2fcTemplate = loadFile(options_intern.path_to_COM2POSE + File.separator +
                options_intern.f_report_resources_basicdata_entry_html);

        log2fcTemplate = log2fcTemplate.replace("{NAME}", "LOG2FC");

        log2fcTemplate =
                log2fcTemplate.replace("{DATA}", getTabularData(transcriptionFactor.name, transcriptionFactor.log2fc));

        template = template.replace("{LOG2FC}", log2fcTemplate);

        template = template.replace("{TPM}", getBasicDataEntry("TPM", transcriptionFactor.tpm));

        template = template.replace("{NORMEX}", getBasicDataEntry("Norm. expression", transcriptionFactor.normex));

        return template;
    }

    private String getBasicData(TranscriptionFactorGroup tfGroup) throws IOException
    {
        if (!tfGroup.realGroup)
        {
            return getBasicData(tfGroup.transcriptionFactors.get(0));
        }

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

    private String getButtonBar(TranscriptionFactorGroup tfGroup) throws IOException
    {
        return getButtonBar(tfGroup.name, tfGroup.hasValidation, tfGroup.hasDistribution, tfGroup.hasRegression);
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

    private boolean generateValidation(TranscriptionFactorGroup tfGroup) throws IOException
    {
        String basicData = getBasicData(tfGroup);

        File templateFile = new File(options_intern.path_to_COM2POSE + File.separator +
                options_intern.f_report_resources_validation_validation_html);

        File d_igv_screenshots =
                new File(options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_igv);

        File d_out_validation =
                new File(options_intern.com2pose_working_directory + File.separator + options_intern.d_out_validation);

        File d_heatmaps = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_heatmap);

        String frame = loadFrame();
        frame = frame.replace("{BODY}", loadFile(templateFile.getAbsolutePath()));

        frame = frame.replace("{TFNAME}", tfGroup.name);

        frame = frame.replace("{TITLE}", tfGroup.name + " - Validation");

        frame = frame.replace("{BASICDATA}", basicData);

        {
            File source = null;

            for (File entry : Objects.requireNonNull(d_heatmaps.listFiles()))
            {
                String entryName = entry.getName();
                if (entryName.equals(tfGroup.name))
                {
                    source = entry;
                    break;
                }
            }
            String id = "validationHeatmaps";
            File target =
                    new File(d_out_validation.getAbsolutePath() + File.separator + tfGroup.name + File.separator + id);
            frame = frame.replace("{VALIDATION_HEATMAP}",
                    (source == null) ? "" : generateTwoLevelImageSelector(id, source, target, false));
        }

        {
            File d_own_tf = new File(
                    d_igv_screenshots.getAbsolutePath() + File.separator + options_intern.folder_out_igv_own_data);
            File source = null;

            for (File entry : Objects.requireNonNull(d_own_tf.listFiles()))
            {
                String entryName = entry.getName();
                entryName = entryName.split("_")[1];
                if (entryName.equals(tfGroup.name))
                {
                    source = entry;
                    break;
                }
            }

            String id = "validationOwnTF";
            File target =
                    new File(d_out_validation.getAbsolutePath() + File.separator + tfGroup.name + File.separator + id);
            frame = frame.replace("{VALIDATION_OWN_TF}", (source == null) ? "" :
                    generateThreeLevelImageSelector(id, source, target, new ArrayList<>(List.of(tfGroup.name)), false,
                            true));
        }

        {
            File d_chip_atlas = new File(d_igv_screenshots.getAbsolutePath() + File.separator +
                    options_intern.folder_out_igv_chip_atlas_data);
            File source = null;

            for (File entry : Objects.requireNonNull(d_chip_atlas.listFiles()))
            {
                String entryName = entry.getName();

                if (entryName.matches("[0-9]+_" + tfGroup.name))
                {
                    source = entry;
                    break;
                }
            }


            String id = "validationChipAtlas";
            File target =
                    new File(d_out_validation.getAbsolutePath() + File.separator + tfGroup.name + File.separator + id);

            frame = frame.replace("{VALIDATION_CHIP_ATLAS_DISABLED}", (source == null) ? "disabled" : "");

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
                                relevantFileName = relevantFileName.replace("_" + tfGroup.name, "");

                                File targetFile = new File(
                                        target.getAbsolutePath() + File.separator + d_hm.getName() + File.separator +
                                                d_group.getName() + File.separator + relevantFileName);

                                copyFile(f_plot, targetFile);
                            }
                        }
                    }
                }
                frame = frame.replace("{VALIDATION_CHIP_ATLAS}",
                        generateThreeLevelImageSelector(id, target, null, true));
            }


        }

        {
            File sourceParentDir = new File(options_intern.com2pose_working_directory + File.separator +
                    options_intern.folder_out_distribution + File.separator +
                    options_intern.folder_out_distribution_logos);

            {
                File sourceDir = new File(sourceParentDir.getAbsolutePath() + File.separator +
                        options_intern.folder_out_distribution_logos_biophysiscal_model);

                File targetFile = new File(
                        options_intern.com2pose_working_directory + File.separator + options_intern.d_out_validation +
                                File.separator + tfGroup.name + File.separator + options_intern.d_out_validation_logos +
                                File.separator + options_intern.f_out_validation_logos_biophysical_png);

                File sourceFile = null;

                for (File directory : Objects.requireNonNull(sourceDir.listFiles()))
                {
                    if (directory.getName().matches("[0-9]+_.*") && directory.isDirectory())
                    {
                        if (directory.getName().split("_")[1].equals(tfGroup.name))
                        {
                            for (File file : Objects.requireNonNull(directory.listFiles()))
                            {
                                if (file.getName().endsWith(".png"))
                                {
                                    sourceFile = file;
                                    break;
                                }
                            }
                            break;
                        }
                    }
                }

                if (sourceFile != null)
                {
                    copyFile(sourceFile, targetFile);
                }
            } // BIOPHYSICAL MODEL

            {
                File sDir = new File(sourceParentDir.getAbsolutePath() + File.separator +
                        options_intern.folder_out_distribution_logos_TF_sequence);

                File targetDir = new File(
                        options_intern.com2pose_working_directory + File.separator + options_intern.d_out_validation +
                                File.separator + tfGroup.name + File.separator +
                                options_intern.d_out_validation_logos_tf_sequence);

                File sourceDir = null;

                for (File directory : Objects.requireNonNull(sDir.listFiles()))
                {
                    if (directory.getName().matches("[0-9]+_.*") && directory.isDirectory())
                    {
                        if (directory.getName().split("_")[1].equals(tfGroup.name))
                        {
                            sourceDir = directory;
                            break;
                        }
                    }
                }

                if (sourceDir != null)
                {
                    ArrayList<File> files = new ArrayList<>();
                    for (File file : Objects.requireNonNull(sourceDir.listFiles()))
                    {
                        if (file.getName().endsWith(".svg"))
                        {
                            File targetFile = new File(targetDir.getAbsolutePath() + File.separator + file.getName());
                            copyFile(file, targetFile);
                            files.add(targetFile);
                        }
                    }

                    if (files.size() > 0)
                    {
                        StringBuilder sb_tf_sequence = new StringBuilder();

                        if (files.size() > 1)
                        {
                            sb_tf_sequence.append("<div class='buttonbar'>");
                            boolean first = true;
                            for (File file : files)
                            {
                                sb_tf_sequence.append(
                                                "<button onclick=\"selectImage(this, 'tf-sequence-image')\" value='" +
                                                        options_intern.d_out_validation_logos_tf_sequence + File.separator +
                                                        file.getName() + "' class='selector").append(first ? " active" : "")
                                        .append("'>");
                                first = false;
                                sb_tf_sequence.append(file.getName(), 0, file.getName().lastIndexOf("."));
                                sb_tf_sequence.append("</button>");
                            }
                            sb_tf_sequence.append("</div>");
                        }

                        sb_tf_sequence.append("<img id='tf-sequence-image' src='" +
                                options_intern.d_out_validation_logos_tf_sequence + File.separator +
                                files.get(0).getName() + "'>");

                        frame = frame.replace("{TF_SEQUENCE}", sb_tf_sequence.toString());
                    }
                }
            } // TF Sequence

            {
                File sDir = new File(sourceParentDir.getAbsolutePath() + File.separator +
                        options_intern.folder_out_distribution_logos_binding_sequence);

                File targetDir = new File(
                        options_intern.com2pose_working_directory + File.separator + options_intern.d_out_validation +
                                File.separator + tfGroup.name + File.separator +
                                options_intern.d_out_validation_logos_tf_binding_sequence);

                File sourceDir = null;

                for (File directory : Objects.requireNonNull(sDir.listFiles()))
                {
                    if (directory.getName().matches("[0-9]+_.*") && directory.isDirectory())
                    {
                        if (directory.getName().split("_")[1].equals(tfGroup.name))
                        {
                            sourceDir = directory;
                            break;
                        }
                    }
                }

                if (sourceDir != null)
                {
                    StringBuilder sb_tf_binding_sequences = new StringBuilder();
                    sb_tf_binding_sequences.append("<div class='buttonbar'>");

                    boolean first = true;
                    String firstImage = "";

                    for (File imageFile : Objects.requireNonNull(sourceDir.listFiles()))
                    {
                        if (!imageFile.getName().endsWith(".png"))
                        {
                            continue;
                        }

                        String relevantFileName = imageFile.getName().split("_")[1];

                        File targetFile =
                                new File(targetDir.getAbsolutePath() + File.separator + relevantFileName + ".png");

                        copyFile(imageFile, targetFile);

                        sb_tf_binding_sequences.append("<button class='selector" + (first ? " active" : "") +
                                "' onclick=\"selectImage(this, 'tf-binding-sequence-image')\" value='" +
                                options_intern.d_out_validation_logos_tf_binding_sequence + File.separator +
                                relevantFileName + ".png'>");
                        sb_tf_binding_sequences.append(relevantFileName);
                        sb_tf_binding_sequences.append("</button>");

                        if (first)
                        {
                            firstImage = relevantFileName;
                            first = false;
                        }
                    }
                    sb_tf_binding_sequences.append("</div>");

                    sb_tf_binding_sequences.append("<img id='tf-binding-sequence-image' src='" +
                            options_intern.d_out_validation_logos_tf_binding_sequence + File.separator + firstImage +
                            ".png'>");

                    frame = frame.replace("{TF_BINDING_SEQUENCE}", sb_tf_binding_sequences.toString());
                }
            } // TF binding sequence
        } // LOGOS

        {
            File sourceDir = new File(
                    options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_igv +
                            File.separator + options_intern.folder_out_igv_dcg_target_genes);

            File targetDir = new File(
                    options_intern.com2pose_working_directory + File.separator + options_intern.d_out_validation);

            File source = null;

            for (File entry : Objects.requireNonNull(sourceDir.listFiles()))
            {
                String entryName = entry.getName();

                if (entryName.equals(tfGroup.name))
                {
                    source = entry;
                    break;
                }
            }

            String id = "validationIGV";
            File target = new File(targetDir.getAbsolutePath() + File.separator + tfGroup.name + File.separator + id);

            frame = frame.replace("{VALIDATION_IGV_DISABLED}", (source == null) ? "disabled" : "");

            frame = frame.replace("{VALIDATION_IGV}", (source == null) ? "" :
                    generateThreeLevelImageSelector(id, source, target, new ArrayList<>(List.of(tfGroup.name)), false,
                            true));
        } // IGV

        frame = relativate(frame, 2);

        {
            StringBuilder sb_actions = new StringBuilder();

            for (TranscriptionFactor tf : tfGroup.transcriptionFactors)
            {
                String command =
                        "window.open('" + options_intern.link_report_genecards.replace("{GENE}", tf.name) + "');\n";
                sb_actions.append(command);
            }

            frame = frame.replace("{GENECARD_BUTTON_ACTION}", sb_actions.toString());
        }

        writeFile(options_intern.com2pose_working_directory + File.separator + options_intern.d_out_validation +
                File.separator + tfGroup.name + File.separator + tfGroup.name + ".html", frame);

        return true;
    }

    private String generateThreeLevelImageSelector(String name, File sourceDir, File targetDir,
                                                   boolean keepFileNameAsIs) throws IOException
    {
        return generateThreeLevelImageSelector(name, sourceDir, targetDir, new ArrayList<>(), keepFileNameAsIs, true);
    }

    private String generateThreeLevelImageSelector(String name, File sourceDir, File targetDir,
                                                   ArrayList<String> specialRemovables, boolean keepFileNameAsIs,
                                                   boolean compression) throws IOException
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
                    ArrayList<String> removables = new ArrayList<>(List.of(group, subgroup, "threshold"));
                    removables.addAll(specialRemovables);

                    if (image_file.isDirectory() || !image_file.getName().endsWith(suffix))
                    {
                        continue;
                    }

                    String relevantFileName = image_file.getName();
                    relevantFileName = relevantFileName.replace(suffix, "");

                    if (!keepFileNameAsIs)
                    {
                        for (String entry : removables)
                        {
                            relevantFileName = relevantFileName.replace(entry, "");
                        }


                        relevantFileName = relevantFileName.replaceAll("_+", "_");

                        while (relevantFileName.endsWith("_"))
                        {
                            relevantFileName = relevantFileName.substring(0, relevantFileName.length() - 1);
                        }

                        while (relevantFileName.startsWith("_"))
                        {
                            relevantFileName = relevantFileName.substring(1);
                        }
                    }

                    if (targetDir != null)
                    {
                        copyFile(image_file, new File(
                                targetDir.getAbsolutePath() + File.separator + group + File.separator + subgroup +
                                        File.separator + relevantFileName + suffix), compression);
                    }

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
                    ArrayList<String> genes = combinations.get(hm).get(group);
                    genes.sort(new StringComparator());

                    for (String gene : genes)
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

    private String generateTwoLevelImageSelector(String name, File sourceDir, File targetDir, boolean compression)
            throws IOException
    {
        String suffix = ".png";
        HashMap<String, ArrayList<String>> combinations = new HashMap<>();

        for (File group_dir : Objects.requireNonNull(sourceDir.listFiles()))
        {
            if (group_dir.isFile())
            {
                continue;
            }

            String group = group_dir.getName();

            combinations.put(group, new ArrayList<>());


            for (File image_file : Objects.requireNonNull(group_dir.listFiles()))
            {
                if (image_file.isDirectory() || !image_file.getName().endsWith(suffix))
                {
                    continue;
                }

                String relevantFileName = image_file.getName();
                relevantFileName = relevantFileName.replace(suffix, "");

                if (targetDir != null)
                {
                    copyFile(image_file, new File(
                            targetDir.getAbsolutePath() + File.separator + group + File.separator + relevantFileName +
                                    suffix), compression);
                }

                combinations.get(group).add(relevantFileName);
            }

        }

        String two_level_image_selector = loadFile(options_intern.f_report_resources_two_level_image_selector_html);

        Set<String> groups = combinations.keySet();
        Set<String> subgroups = new HashSet<>();

        {
            StringBuilder sb_groups = new StringBuilder();

            for (String group : groups)
            {
                subgroups.addAll(combinations.get(group));

                sb_groups.append("<button class=\"{ID} group-selector\" onclick=\"select_group" + "('{ID}', this, " +
                        "{ID}Combinations)\" " + "value=\"" + group + "\">" + group + "</button>");
            }
            two_level_image_selector = two_level_image_selector.replace("{GROUPS}", sb_groups.toString());
        }

        {
            StringBuilder sb_subgroups = new StringBuilder();
            for (String subgroup : subgroups)
            {
                sb_subgroups.append(
                        "<button class=\"{ID} subgroup-selector\" onclick=\"select_subgroup" + "('{ID}', " + "this, " +
                                "{ID}Combinations)\" " + "value=\"" + subgroup + "\">" + subgroup + "</button>");
            }

            two_level_image_selector = two_level_image_selector.replace("{SUBGROUPS}", sb_subgroups.toString());
        }

        String json;
        {
            HashMap<String, String> lv1 = new HashMap<>();

            for (String group : combinations.keySet())
            {
                StringBuilder sb_subgroups = new StringBuilder("[");
                for (String subgroup : combinations.get(group))
                {
                    sb_subgroups.append("\"");
                    sb_subgroups.append(subgroup);
                    sb_subgroups.append("\",");
                }
                sb_subgroups.setLength(sb_subgroups.length() - 1);
                sb_subgroups.append("]");
                lv1.put(group, sb_subgroups.toString());
            }
            json = mapToJson(lv1);
        }

        two_level_image_selector = two_level_image_selector.replace("{ID}", name);

        two_level_image_selector = two_level_image_selector.replace("{COMBINATIONS}", json);

        return two_level_image_selector;
    }

    private boolean generateDistribution(TranscriptionFactorGroup tfGroup) throws IOException
    {
        String basicData = getBasicData(tfGroup);

        File templateFile = new File(options_intern.path_to_COM2POSE + File.separator +
                options_intern.f_report_resources_distribution_distribution_html);

        File d_distribution_output = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.d_out_distribution +
                        File.separator + tfGroup.name);

        ArrayList<String> existingHMs = new ArrayList<>();

        File d_plots_hm = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_distribution +
                        File.separator + options_intern.folder_out_distribution_plots + File.separator +
                        options_intern.folder_out_distribution_plots_HM);

        for (File d_hm : Objects.requireNonNull(d_plots_hm.listFiles()))
        {
            for (File f_plot : Objects.requireNonNull(d_hm.listFiles()))
            {
                if (f_plot.getName().substring(0, f_plot.getName().lastIndexOf(".")).equals(tfGroup.name))
                {
                    existingHMs.add(d_hm.getName());
                    copyFile(f_plot, new File(
                            d_distribution_output.getAbsolutePath() + File.separator + d_hm.getName() + ".png"));
                    break;
                }
            }
        }

        File d_plots_all = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_distribution +
                        File.separator + options_intern.folder_out_distribution_plots + File.separator +
                        options_intern.folder_out_distribution_plots_ALL);

        for (File f_plot : Objects.requireNonNull(d_plots_all.listFiles()))
        {
            if (f_plot.getName().substring(0, f_plot.getName().lastIndexOf(".")).equals(tfGroup.name))
            {
                String categoryName = d_plots_all.getName().substring(d_plots_all.getName().indexOf("_") + 1);
                existingHMs.add(categoryName);
                copyFile(f_plot,
                        new File(d_distribution_output.getAbsolutePath() + File.separator + categoryName + ".png"));
                break;
            }
        }

        if (existingHMs.isEmpty())
        {
            return false;
        }

        String frame = loadFrame();
        frame = frame.replace("{BODY}", loadFile(templateFile.getAbsolutePath()));

        frame = frame.replace("{TFNAME}", tfGroup.name);
        frame = frame.replace("{BASICDATA}", basicData);

        frame = frame.replace("{TITLE}", tfGroup.name + " - Distribution");

        StringBuilder sb_histoneModifications = new StringBuilder();

        boolean first = true;

        for (String histoneModification : existingHMs)
        {
            sb_histoneModifications.append(
                    "<button onclick=\"selectImage(this, 'distribution-plot')\" " + "class=\"selector" +
                            (first ? " active" : "") + "\" value=\"" + histoneModification + ".png\">" +
                            histoneModification + "</button>");
            first = false;
        }

        frame = frame.replace("{HISTONEMODIFICATIONS}", sb_histoneModifications.toString());
        frame = frame.replace("{DISTRIBUTION-PLOT-PATH}", existingHMs.get(0) + ".png");

        {
            StringBuilder sb_actions = new StringBuilder();

            for (TranscriptionFactor tf : tfGroup.transcriptionFactors)
            {
                String command =
                        "window.open('" + options_intern.link_report_genecards.replace("{GENE}", tf.name) + "');\n";
                sb_actions.append(command);
            }

            frame = frame.replace("{GENECARD_BUTTON_ACTION}", sb_actions.toString());
        }

        {
            File allFile = new File(options_intern.com2pose_working_directory + File.separator +
                    options_intern.folder_out_distribution + File.separator +
                    options_intern.folder_out_distribution_stats + File.separator +
                    options_intern.folder_out_distribution_stats_ALL + File.separator +
                    options_intern.file_suffix_distribution_analysis_plot_stats);
            File hmsDir = new File(options_intern.com2pose_working_directory + File.separator +
                    options_intern.folder_out_distribution + File.separator +
                    options_intern.folder_out_distribution_stats + File.separator +
                    options_intern.folder_out_distribution_stats_HM);
            HashMap<String, Integer> ranks = new HashMap<>();
            HashMap<String, Integer> sizes = new HashMap<>();

            try
            {
                ranks.put(options_intern.distribution_analysis_all_name,
                        Integer.parseInt(findValueInTable(tfGroup.name, 1, 0, allFile, "\t", true)));
            } catch (NoSuchFieldException e)
            {
                ranks.put(options_intern.distribution_analysis_all_name, -1);
            } finally
            {
                String[] lines = loadFile(allFile.getAbsolutePath()).split("\n");
                String lastLine = lines[lines.length - 1];
                sizes.put(options_intern.distribution_analysis_all_name, Integer.parseInt(lastLine.split("\t")[0]));
            }

            for (File hmDir : Objects.requireNonNull(hmsDir.listFiles()))
            {
                if (!hmDir.isDirectory())
                {
                    continue;
                }

                File statsFile = new File(hmDir.getAbsolutePath() + File.separator +
                        options_intern.file_suffix_distribution_analysis_plot_stats);

                if (statsFile.exists())
                {
                    try
                    {
                        ranks.put(hmDir.getName(),
                                Integer.parseInt(findValueInTable(tfGroup.name, 1, 0, statsFile, "\t", true)));
                    } catch (NoSuchFieldException e)
                    {
                        ranks.put(hmDir.getName(), -1);
                    } finally
                    {
                        String[] lines = loadFile(statsFile.getAbsolutePath()).split("\n");
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
                sb_dcg.append(
                        entry.getValue() == -1 ? "{RELATIVATION}not_available.png" : "{RELATIVATION}is_available.png");
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

            frame = frame.replace("{DISCOUNTED_CUMULATIVE_GAIN}", sb_dcg.toString());
        }

        frame = relativate(frame, 2);

        writeFile(d_distribution_output + File.separator + tfGroup.name + ".html", frame);

        return true;
    }


    private boolean generateRegression(TranscriptionFactorGroup tfGroup) throws IOException
    {
        File templateFile = new File(options_intern.path_to_COM2POSE + File.separator +
                options_intern.f_report_resources_regression_regression_html);

        File d_in_plots =
                new File(options_intern.com2pose_working_directory + File.separator + options_intern.folder_plots);

        File d_out_regression =
                new File(options_intern.com2pose_working_directory + File.separator + options_intern.d_out_regression);

        String frame = loadFrame();

        frame = frame.replace("{BODY}", loadFile(templateFile.getAbsolutePath()));

        frame = frame.replace("{TFNAME}", tfGroup.name);

        frame = frame.replace("{TITLE}", tfGroup.name + " - Regression");

        File tfDir = new File(d_out_regression.getAbsolutePath() + File.separator + tfGroup.name);

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

                    copyFile(image_file, targetFile);
                }
            }
        }

        String overviewCoefficients = generateThreeLevelImageSelector(overviewCoefficientsID,
                new File(tfDir.getAbsolutePath() + File.separator + overviewCoefficientsID), null, true);

        String overviewHeatmaps = generateThreeLevelImageSelector(overviewHeatmapsID,
                new File(tfDir.getAbsolutePath() + File.separator + overviewHeatmapsID), null, true);

        frame = frame.replace("{COEFFICIENTS}",
                getTabularData("regressionCoefficients", tfGroup.regressionCoefficients));

        frame = frame.replace("{OVERVIEW_COEFFICIENTS}", overviewCoefficients);

        frame = frame.replace("{OVERVIEW_HEATMAPS}", overviewHeatmaps);

        frame = frame.replace("{ADDHEAD}", "<script src=\"COMBINATIONS.js\"></script>");

        frame = relativate(frame, 2);

        writeFile(options_intern.com2pose_working_directory + File.separator + options_intern.d_out_regression +
                File.separator + tfGroup.name + File.separator + tfGroup.name + ".html", frame);

        return true;
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

        copyFile(
                new File(options_intern.path_to_COM2POSE + File.separator + options_intern.f_report_resources_logo_png),
                new File(options_intern.com2pose_working_directory + File.separator +
                        options_intern.f_out_report_logo_png));

        copyFile(new File(
                        options_intern.path_to_COM2POSE + File.separator + options_intern.f_report_resources_is_available_png),
                new File(options_intern.com2pose_working_directory + File.separator +
                        options_intern.f_out_report_is_available_png));

        copyFile(new File(
                        options_intern.path_to_COM2POSE + File.separator + options_intern.f_report_resources_not_available_png),
                new File(options_intern.com2pose_working_directory + File.separator +
                        options_intern.f_out_report_not_available_png));
    }

    private void copyFile(File source, File target) throws IOException
    {
        copyFile(source, target, true);
    }

    private void copyFile(File source, File target, boolean compression) throws IOException
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

            if (target.getName().endsWith(".png") && compression)
            {
                String command = "pngtopnm " + source.getAbsolutePath() + " | pnmquant 16 | pnmtopng > " +
                        target.getAbsolutePath();
                String[] cmd = {"/bin/sh", "-c", command};
                executorService.submit(() ->
                {
                    try
                    {
                        Process child = Runtime.getRuntime().exec(cmd);
                        child.waitFor();

                    } catch (IOException | InterruptedException e)
                    {
                        e.printStackTrace();
                    }

                });
            } else
            {
                executorService.submit(() ->
                {
                    try
                    {
                        Files.copy(source.toPath(), target.toPath(), REPLACE_EXISTING);
                    } catch (IOException e)
                    {
                        e.printStackTrace();
                    }
                });
            }
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

    private void generateImportantLoci() throws IOException
    {
        File sourceDir = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_igv +
                        File.separator + options_intern.folder_out_igv_important_loci);

        File targetDir = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.d_out_important_loci);


        for (File groupFile : Objects.requireNonNull(sourceDir.listFiles()))
        {
            for (File imageFile : Objects.requireNonNull(groupFile.listFiles()))
            {
                if (!imageFile.getName().endsWith(".png"))
                {
                    continue;
                }

                for (String importantTf : options_intern.igv_important_locus_all_prio_tf)
                {
                    if (imageFile.getName().startsWith(importantTf))
                    {
                        File targetFile = new File(
                                targetDir.getAbsolutePath() + File.separator + groupFile.getName() + File.separator +
                                        importantTf + File.separator + imageFile.getName());

                        copyFile(imageFile, targetFile);
                        break;
                    }
                }
            }
        }

        String frame = loadFrame();
        String important_loci = frame.replace("{BODY}", loadFile(options_intern.path_to_COM2POSE + File.separator +
                options_intern.f_report_resources_important_loci_html));

        important_loci = important_loci.replace("{TITLE}", "Important loci");
        important_loci = important_loci.replace("{IMAGES}",
                generateThreeLevelImageSelector(targetDir.getName(), targetDir, null, true));
        important_loci = relativate(important_loci, 0);

        writeFile(options_intern.com2pose_working_directory + File.separator +
                options_intern.f_out_report_important_loci_html, important_loci);
    }

    private void generateTopLog2fc() throws IOException
    {
        File sourceDir = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_igv +
                        File.separator + options_intern.folder_out_igv_top_log2fc);

        File targetDir =
                new File(options_intern.com2pose_working_directory + File.separator + options_intern.d_out_top_log2fc);

        String frame = loadFrame();
        String top_log2fc = frame.replace("{BODY}", loadFile(
                options_intern.path_to_COM2POSE + File.separator + options_intern.f_report_resources_top_log2fc_html));

        top_log2fc = top_log2fc.replace("{TITLE}", "Top log2fc");
        top_log2fc = top_log2fc.replace("{IMAGES}",
                generateThreeLevelImageSelector(targetDir.getName(), sourceDir, targetDir, false));
        top_log2fc = relativate(top_log2fc, 0);

        writeFile(options_intern.com2pose_working_directory + File.separator +
                options_intern.f_out_report_top_log2fc_html, top_log2fc);
    }

    private void generateCoOccurrence() throws IOException
    {
        File dataSource = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_distribution +
                        File.separator + options_intern.folder_out_distribution_cooccurence + File.separator +
                        options_intern.file_suffix_cooccurence_frequencies);

        String input = loadFile(dataSource.getAbsolutePath());
        String frame = loadFrame();

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

        String cooccurrence = frame.replace("{BODY}", loadFile(options_intern.path_to_COM2POSE + File.separator +
                options_intern.f_report_resources_cooccurrence_html));

        cooccurrence = cooccurrence.replace("{TABLE}", getTabularData("coOccurrence", data));
        cooccurrence = cooccurrence.replace("{TITLE}", "Co-Occurrence analysis");

        cooccurrence = relativate(cooccurrence, 0);

        writeFile(options_intern.com2pose_working_directory + File.separator +
                options_intern.f_out_report_cooccurrence_html, cooccurrence);
    }

    static class StringComparator implements Comparator<String>
    {
        @Override public int compare(String a, String b)
        {
            return prefixNum(a) - prefixNum(b);
        }

        private int prefixNum(String a)
        {
            if (a.matches("[0-9]+_.*"))
            {
                return Integer.parseInt(a.split("_")[0]);
            } else
            {
                return 0;
            }
        }
    }
}
