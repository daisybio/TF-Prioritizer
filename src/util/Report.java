package util;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.*;

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

    private record TranscriptionFactor(String geneID, String name, Map<String, Number> log2fc, Map<String, Number> tpm,
                                       Map<String, Number> normex, ArrayList<String> histoneModifications,
                                       boolean hasIGV)
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
                        Map<String, Number> log2fc = new HashMap<>();
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
                                        String name = group1 + " | " + group2;
                                        log2fc.put(name, log2fc_value);
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

        for (Map.Entry<String, Number> kvPair : data.entrySet())
        {
            sb_data.append("<div class=\"keyvaluepair\"><p>" + kvPair.getKey() + "</p><p>" +
                    formatter.format(kvPair.getValue()) + "</p></div>");
        }

        template = template.replace("{DATA}", sb_data.toString());

        return template;
    }

    private String getBasicData(TranscriptionFactor transcriptionFactor) throws IOException
    {
        String template = loadFile(
                options_intern.path_to_COM2POSE + File.separator + options_intern.f_report_resources_basicdata_html);


        template = template.replace("{LOG2FC}", getBasicDataEntry("LOG2FC", transcriptionFactor.log2fc));

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
        return getButtonBar(tf.name, tf.hasIGV, true, false);
    }

    private String getButtonBar(TranscriptionFactorGroup tfGroup) throws IOException
    {
        return getButtonBar(tfGroup.name, tfGroup.hasIGV, true, false);
    }

    private String getButtonBar(String name, boolean hasValidation, boolean hasDistribution, boolean hasRegression)
            throws IOException
    {
        String buttonbar = loadFile(options_intern.path_to_COM2POSE + File.separator +
                options_intern.f_report_resources_home_buttonbar_html);

        buttonbar = buttonbar.replace("{VALIDATION}", "VALIDATION/" + name + ".html");
        buttonbar = buttonbar.replace("{DISTRIBUTION}", "DISTRIBUTION/" + name + ".html");
        buttonbar = buttonbar.replace("{REGRESSION}", "REGRESSION/" + name + ".html");

        buttonbar = buttonbar.replace("{HASVALIDATION}", hasValidation ? "" : "disabled");
        buttonbar = buttonbar.replace("{HASREGRESSION}", hasRegression ? "" : "disabled");
        buttonbar = buttonbar.replace("{HASDISTRIBUTION}", hasDistribution ? "" : "disabled");

        return buttonbar;
    }


    private void generateValidation(TranscriptionFactorGroup tfGroup) throws IOException
    {
        if (tfGroup.realGroup)
        {
            generateValidation(tfGroup.name, getBasicData(tfGroup), new ArrayList<String>(),
                    new HashMap<String, Number>());
        } else
        {
            TranscriptionFactor tf = tfGroup.transcriptionFactors.get(0);

            generateValidation(tf.name, getBasicData(tf), tf.histoneModifications, tf.log2fc);
        }
    }

    private void generateValidation(String name, String basicData, ArrayList<String> histoneModifications,
                                    Map<String, Number> log2fc) throws IOException
    {
        File templateFile = new File(options_intern.path_to_COM2POSE + File.separator +
                options_intern.f_report_resources_validation_validation_html);

        String templateFrame = loadFrame();
        templateFrame = templateFrame.replace("{BODY}", loadFile(templateFile.getAbsolutePath()));


        String frame = templateFrame;

        frame = frame.replace("{TFNAME}", name);

        frame = frame.replace("{TITLE}", name + " - Validation");

        frame = frame.replace("{BASICDATA}", basicData);

        StringBuilder sb_histoneModifications = new StringBuilder();

        {   // Histone modifications
            for (String hisoneModification : histoneModifications)
            {
                sb_histoneModifications.append(
                        "<a href=\"{RELATIVATION}PARAMETERS.html\" " + "class=\"button\">" + hisoneModification +
                                "</a>");
            }

            frame = frame.replace("{HISTONEMODIFICATIONS}", sb_histoneModifications.toString());
        }   // Histone modifications

        {   // Groups
            StringBuilder sb_groups = new StringBuilder();

            for (Map.Entry<String, Number> group : log2fc.entrySet())
            {
                sb_groups.append(
                        "<a href=\"{RELATIVATION}PARAMETERS.html\" " + "class=\"button\">" + group.getKey() + "</a>");
            }

            frame = frame.replace("{GROUPS}", sb_groups.toString());
        }   // Groups

        frame = relativate(frame, 1);

        writeFile(options_intern.com2pose_working_directory + File.separator + options_intern.d_out_validation +
                File.separator + name + ".html", frame);

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

        String templateFrame = loadFrame();
        templateFrame = templateFrame.replace("{BODY}", loadFile(templateFile.getAbsolutePath()));


        String frame = templateFrame;

        frame = frame.replace("{TFNAME}", name);
        frame = frame.replace("{BASICDATA}", basicData);

        frame = frame.replace("{TITLE}", name + " - Distribution");

        frame = relativate(frame, 1);

        writeFile(options_intern.com2pose_working_directory + File.separator + options_intern.d_out_distribution +
                File.separator + name + ".html", frame);
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

        String templateFrame = loadFrame();
        templateFrame = templateFrame.replace("{BODY}", loadFile(templateFile.getAbsolutePath()));

        String frame = templateFrame;

        frame = frame.replace("{TFNAME}", name);

        frame = frame.replace("{TITLE}", name + " - Regression");

        frame = relativate(frame, 1);

        writeFile(options_intern.com2pose_working_directory + File.separator + options_intern.d_out_regression +
                File.separator + name + ".html", frame);

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

    private void writeFile(String path, String content) throws IOException
    {
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
}
