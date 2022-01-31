package util.Report;

import util.Logger;
import util.Options_intern;

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
    static Logger logger = null;
    static Options_intern options_intern = null;
    final ArrayList<TranscriptionFactorGroup> transcriptionFactorGroups = new ArrayList<>();
    static final DecimalFormat formatter = new DecimalFormat("0.###");
    static Map<SelectorTypes, ArrayList<String>> existingValues = new HashMap<>();

    public Report(Options_intern options_intern) throws IOException
    {
        Report.options_intern = options_intern;

        for (SelectorTypes type : SelectorTypes.values())
        {
            existingValues.put(type, new ArrayList<>());
        }

        logger = new Logger(true, options_intern.com2pose_working_directory);

        logger.logLine("[REPORT] Start loading TF data");
        loadTFs();
        findExistingValues();
        logger.logLine("[REPORT] Finished loading TF data");
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

                String data = FileManagement.loadFile(f_data.getAbsolutePath());

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

        existingValues.get(SelectorTypes.HISTONE_MODIFICATIONS).addAll(histoneModifications);

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
                        String geneID = FileManagement.findValueInTable(tf_name, 1, 0, geneIDFile, "\t", true);
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

                                    existingValues.get(SelectorTypes.GROUPS).add(group1);
                                    existingValues.get(SelectorTypes.GROUPS).add(group2);

                                    try
                                    {
                                        double log2fc_value = Double.parseDouble(
                                                FileManagement.findValueInTable(geneID, 0, 1, entry, "\t", false));

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
                                        double tpm_value = Double.parseDouble(
                                                FileManagement.findValueInTable(geneID, 0, 3, entry, "\t", false));
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
                                        int exp_value = Integer.parseInt(
                                                FileManagement.findValueInTable(geneID, 1, 2, entry, "\t", false));
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
                            new TranscriptionFactorGroup(tfGroupName, tf_group, regressionCoefficients));
                }
            }
        }
    }

    public void generate() throws IOException, InterruptedException
    {
        logger.logLine("[REPORT] Start generating report");

        copyDependencies();

        PageGenerators.generateGeneralPages();
        PageGenerators.generateHome(transcriptionFactorGroups);

        int i = 1;

        for (TranscriptionFactorGroup tfGroup : transcriptionFactorGroups)
        {
            System.out.print("Generating report for tf " + i + " out of " + transcriptionFactorGroups.size() + ": " +
                    tfGroup.getName() + ": Validation.\r");
            tfGroup.setValidation(PageGenerators.generateValidation(tfGroup));
            System.out.print("Generating report for tf " + i + " out of " + transcriptionFactorGroups.size() + ": " +
                    tfGroup.getName() + ": Distribution.\r");
            tfGroup.setDistribution(PageGenerators.generateDistribution(tfGroup));
            System.out.print("Generating report for tf " + i + " out of " + transcriptionFactorGroups.size() + ": " +
                    tfGroup.getName() + ": Regression.\r");
            tfGroup.setRegression(PageGenerators.generateRegression(tfGroup));
            i++;
        }

        logger.logLine("[REPORT] Finished generating report");
    }

    private void copyDependencies() throws IOException
    {
        logger.logLine("[REPORT] Start copying dependencies");
        String css = FileManagement.loadFile(
                options_intern.path_to_COM2POSE + File.separator + options_intern.f_report_resources_style);
        String script = FileManagement.loadFile(
                options_intern.path_to_COM2POSE + File.separator + options_intern.f_report_resources_script);

        FileManagement.writeFile(
                options_intern.com2pose_working_directory + File.separator + options_intern.f_out_report_style, css);
        FileManagement.writeFile(
                options_intern.com2pose_working_directory + File.separator + options_intern.f_out_report_script,
                script);

        FileManagement.copyFile(
                new File(options_intern.path_to_COM2POSE + File.separator + options_intern.f_report_resources_logo_png),
                new File(options_intern.com2pose_working_directory + File.separator +
                        options_intern.f_out_report_logo_png));

        FileManagement.copyDirectory(
                new File(options_intern.path_to_COM2POSE + File.separator + options_intern.d_report_resources_media),
                new File(options_intern.com2pose_working_directory + File.separator + options_intern.d_out_media),
                false);
        logger.logLine("[REPORT] Finished copying dependencies");
    }

    private void findExistingValues()
    {
        {
            File pairings = new File(options_intern.com2pose_working_directory + File.separator +
                    options_intern.folder_name_deseq2_preprocessing + File.separator +
                    options_intern.folder_name_deseq2_preprocessing_combined);

            for (File pairing : Objects.requireNonNull(pairings.listFiles()))
            {
                existingValues.get(SelectorTypes.GROUP_PAIRINGS).add(pairing.getName());
            }
        }

        existingValues.get(SelectorTypes.IMPORTANT_LOCI).addAll(options_intern.igv_important_locus_all_prio_tf);

        existingValues.get(SelectorTypes.TOP_LOG2FC).addAll(Arrays.asList("downregulated", "upregulated"));
        existingValues.get(SelectorTypes.REGRESSION_CUTOFFS)
                .addAll(Arrays.asList("0.1", "0.2", "0.3", "0.4", "0.5", "0.6"));

        existingValues.get(SelectorTypes.DISTRIBUTION_OPTIONS)
                .addAll(existingValues.get(SelectorTypes.HISTONE_MODIFICATIONS));
        existingValues.get(SelectorTypes.DISTRIBUTION_OPTIONS).add("ALL");

        for (SelectorTypes type : SelectorTypes.values())
        {
            if (existingValues.get(type).size() > 0)
            {
                Set<String> set = new HashSet<>(existingValues.get(type));
                existingValues.get(type).clear();
                existingValues.get(type).addAll(set);
                Collections.sort(existingValues.get(type));
            }
        }
    }
}
