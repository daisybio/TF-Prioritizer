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
    private final ArrayList<Map<String, String>> transcriptionFactors = new ArrayList<>();

    public Report(Options_intern options_intern) throws IOException
    {
        this.options_intern = options_intern;
        logger = new Logger(true, options_intern.com2pose_working_directory);
        loadTFs();
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

        try (Scanner scanner = new Scanner(tf_file))
        {
            boolean firstLine = true;
            while (scanner.hasNextLine())
            {
                String line = scanner.nextLine();
                if (firstLine)
                {
                    firstLine = false;
                    continue;
                }
                String tf_name = line.split("\t")[1];
                try
                {
                    String geneID = findValueInTable(tf_name, 1, 0, geneIDFile, "\t", true);
                    Map<String, String> hashMap = new HashMap<>();
                    hashMap.put("ENSG", geneID.toUpperCase());
                    hashMap.put("GeneSymbol", tf_name.toUpperCase());
                    transcriptionFactors.add(hashMap);
                } catch (NoSuchFieldException ignore)
                {
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
        generateValidation();
        generateDistribution();
        generateRegression();

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
        DecimalFormat formatter = new DecimalFormat("0.000");


        int i = 1;

        for (Map<String, String> tfMap : transcriptionFactors)
        {
            String tf_string = loadFile(
                    options_intern.path_to_COM2POSE + File.separator + options_intern.f_report_resources_home_tf_html);
            String gene_symbol = tfMap.get("GeneSymbol").toUpperCase();
            String ensg_symbol = tfMap.get("ENSG");

            tf_string = tf_string.replace("{TF_NAME}", i + ". " + gene_symbol);

            {   //LOG2FC
                File d_log2fs = new File(options_intern.com2pose_working_directory + File.separator +
                        options_intern.folder_name_deseq2_output);

                StringBuilder sb_log2fc = new StringBuilder();

                for (File entry : Objects.requireNonNull(d_log2fs.listFiles()))
                {
                    if (entry.isFile())
                    {
                        String group1 = entry.getName().split("_")[0];
                        String group2 = entry.getName().split("_")[1];

                        try
                        {
                            double log2fc = Double.parseDouble(findValueInTable(ensg_symbol, 0, 1, entry, "\t", false));
                            sb_log2fc.append("<div class=\"keyvaluepair\"><p>" + group1 + " -> " + group2 + ":</p><p>" +
                                    formatter.format(log2fc) + "</p></div>");
                        } catch (NoSuchFieldException ignored)
                        {
                        }
                    }
                }
                tf_string = tf_string.replace("{LOG2FC}", sb_log2fc.toString());
            }   //LOG2FC

            {   //TPM
                File d_tpm = new File(options_intern.com2pose_working_directory + File.separator +
                        options_intern.folder_name_deseq2_preprocessing + File.separator +
                        options_intern.folder_name_deseq2_preprocessing_tpm + File.separator +
                        options_intern.folder_name_deseq2_preprocessing_tpm_results);

                StringBuilder sb_tpm = new StringBuilder();

                for (File entry : Objects.requireNonNull(d_tpm.listFiles()))
                {
                    if (entry.isFile() && entry.getName().endsWith(".csv"))
                    {
                        String group = entry.getName().split("_")[0];

                        try
                        {
                            double tpm = Double.parseDouble(findValueInTable(ensg_symbol, 0, 3, entry, "\t", false));
                            sb_tpm.append(
                                    "<div class=\"keyvaluepair\"><p>" + group + "</p><p>" + formatter.format(tpm) +
                                            "</p></div>");
                        } catch (NoSuchFieldException ignored)
                        {
                        }
                    }
                }
                tf_string = tf_string.replace("{TPM}", sb_tpm.toString());
            }   //TPM

            {   //Normalized expression
                File d_normex = new File(options_intern.com2pose_working_directory + File.separator +
                        options_intern.folder_name_deseq2_preprocessing + File.separator +
                        options_intern.folder_name_deseq2_preprocessing_gene_symbols);

                StringBuilder sb_normex = new StringBuilder();

                for (File entry : Objects.requireNonNull(d_normex.listFiles()))
                {
                    if (entry.isFile() && entry.getName().endsWith(".csv"))
                    {
                        String group = entry.getName().split("\\.")[0];

                        try
                        {
                            int exp = Integer.parseInt(findValueInTable(ensg_symbol, 1, 2, entry, "\t", false));
                            sb_normex.append(
                                    "<div class=\"keyvaluepair\"><p>" + group + "</p><p>" + exp + "</p" + "></div>");
                        } catch (NoSuchFieldException ignored)
                        {
                        }
                    }
                }
                tf_string = tf_string.replace("{NORMEX}", sb_normex.toString());
            }   //Normalized expression

            tf_string = tf_string.replace("{GENEID}", ensg_symbol);

            tf_string = tf_string.replace("{VALIDATION}", "VALIDATION/" + gene_symbol + ".html");
            tf_string = tf_string.replace("{DISTRIBUTION}", "DISTRIBUTION/" + gene_symbol + ".html");
            tf_string = tf_string.replace("{REGRESSION}", "REGRESSION/" + gene_symbol + ".html");

            sb_tfs.append(tf_string);
            i++;
        }
        home = home.replace("{TFS}", sb_tfs.toString());

        home = relativate(home, 0);


        writeFile(options_intern.com2pose_working_directory + File.separator + options_intern.f_out_report_home, home);
        logger.logLine("[REPORT] Finished generating report home page");
    }

    private void generateValidation() throws IOException
    {
        File templateFile = new File(options_intern.path_to_COM2POSE + File.separator +
                options_intern.f_report_resources_validation_validation_html);


        for (Map<String, String> entry : transcriptionFactors)
        {
            String frame = loadFrame();
            String validation = loadFile(templateFile.getAbsolutePath());
            String tfName = entry.get("GeneSymbol");

            validation = validation.replace("{TFNAME}", tfName);

            frame = frame.replace("{BODY}", validation);

            frame = frame.replace("{TITLE}", tfName + " - Validation");

            frame = relativate(frame, 1);

            writeFile(options_intern.com2pose_working_directory + File.separator + options_intern.d_out_validation +
                    File.separator + tfName + ".html", frame);
        }
    }

    private void generateDistribution() throws IOException
    {
        File templateFile = new File(options_intern.path_to_COM2POSE + File.separator +
                options_intern.f_report_resources_distribution_distribution_html);


        for (Map<String, String> entry : transcriptionFactors)
        {
            String frame = loadFrame();
            String distribution = loadFile(templateFile.getAbsolutePath());
            String tfName = entry.get("GeneSymbol");

            distribution = distribution.replace("{TFNAME}", tfName);

            frame = frame.replace("{BODY}", distribution);

            frame = frame.replace("{TITLE}", tfName + " - Distribution");

            frame = relativate(frame, 1);

            writeFile(options_intern.com2pose_working_directory + File.separator + options_intern.d_out_distribution +
                    File.separator + tfName + ".html", frame);
        }
    }

    private void generateRegression() throws IOException
    {
        File templateFile = new File(options_intern.path_to_COM2POSE + File.separator +
                options_intern.f_report_resources_regression_regression_html);


        for (Map<String, String> entry : transcriptionFactors)
        {
            String frame = loadFrame();
            String regression = loadFile(templateFile.getAbsolutePath());
            String tfName = entry.get("GeneSymbol");

            regression = regression.replace("{TFNAME}", tfName);

            frame = frame.replace("{BODY}", regression);

            frame = frame.replace("{TITLE}", tfName + " - Regression");

            frame = relativate(frame, 1);

            writeFile(options_intern.com2pose_working_directory + File.separator + options_intern.d_out_regression +
                    File.separator + tfName + ".html", frame);
        }
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
