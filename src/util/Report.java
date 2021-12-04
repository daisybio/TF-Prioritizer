package util;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Scanner;

public class Report
{
    private final Logger logger;
    private final Options_intern options_intern;

    public Report(Options_intern options_intern) throws Exception
    {
        this.options_intern = options_intern;
        logger = new Logger(true, options_intern.com2pose_working_directory);
    }

    public void generate() throws Exception
    {
        logger.logLine("[REPORT] Start generating report");

        generateHome();
        styleAndScript();
        generateParameters();

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


        writeFile(options_intern.com2pose_working_directory + File.separator + options_intern.f_out_report_parameters,
                parameters);
        logger.logLine("[REPORT] Finished generating report parameters page");
    }

    private void generateHome() throws Exception
    {
        logger.logLine("[REPORT] Start generating report home page");
        String home = loadFrame();
        home = home.replace("{BODY}", loadFile(
                options_intern.path_to_COM2POSE + File.separator + options_intern.f_report_resources_home_home_html));
        home = home.replace("{TITLE}", "Overview");

        File tf_file = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_distribution +
                        File.separator + options_intern.folder_out_distribution_dcg + File.separator +
                        options_intern.file_suffix_distribution_analysis_dcg);

        try (Scanner scanner = new Scanner(tf_file))
        {
            StringBuilder sb_tfs = new StringBuilder();
            boolean firstLine = true;

            File gene_symbol_map = new File(options_intern.com2pose_working_directory + File.separator +
                    options_intern.folder_name_deseq2_preprocessing + File.separator +
                    options_intern.file_suffix_deseq2_mapping);

            int i = 1;

            while (scanner.hasNextLine())
            {
                String line = scanner.nextLine();
                if (firstLine)
                {
                    firstLine = false;
                    continue;
                }
                String tf_string = loadFile(options_intern.path_to_COM2POSE + File.separator +
                        options_intern.f_report_resources_home_tf_html);
                String gene_symbol = line.split("\t")[1];
                String ensg_symbol = "";

                try (Scanner mapScanner = new Scanner(gene_symbol_map))
                {
                    while (mapScanner.hasNextLine())
                    {
                        String translation = mapScanner.nextLine();
                        if (translation.split("\t").length == 2)
                        {
                            if (gene_symbol.equalsIgnoreCase(translation.split("\t")[1]))
                            {
                                ensg_symbol = translation.split("\t")[0];
                            }
                        }
                    }
                }
                tf_string = tf_string.replace("{TF_NAME}", i + ". " + gene_symbol);
                tf_string = tf_string.replace("{CONTENT}", ensg_symbol);
                sb_tfs.append(tf_string);
                i++;
            }
            home = home.replace("{TFS}", sb_tfs.toString());
        }


        writeFile(options_intern.com2pose_working_directory + File.separator + options_intern.f_out_report_home, home);
        logger.logLine("[REPORT] Finished generating report home page");
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
