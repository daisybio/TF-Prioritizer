package lib.Report;

import tfprio.TFPRIO;
import util.FileManagement;
import util.Logger;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;

public class Report
{
    static Logger logger = null;
    final ArrayList<TranscriptionFactorGroup> transcriptionFactorGroups = new ArrayList<>();
    static final DecimalFormat formatter = new DecimalFormat("0.###");
    static Map<SelectorTypes, ArrayList<String>> existingValues = new HashMap<>();

    public Report() throws IOException
    {
        for (SelectorTypes type : SelectorTypes.values())
        {
            existingValues.put(type, new ArrayList<>());
        }

        logger = new Logger("Report", true, TFPRIO.configs.general.logFile.get());

        logger.info("Start loading TF data");
        loadTFs();
        findExistingValues();
        logger.info("Finished loading TF data");
    }


    private Map<String, Map<String, Map<String, Double>>> loadRegressionCoefficients() throws IOException
    {
        Map<String, Map<String, Map<String, Double>>> coefficients = new HashMap<>();

        File parentDir = TFPRIO.configs.dynamite.fileStructure.d_output.get();

        for (File hm_dir : Objects.requireNonNull(parentDir.listFiles()))
        {
            coefficients.put(hm_dir.getName(), new HashMap<>());

            for (File groups_dir : Objects.requireNonNull(hm_dir.listFiles()))
            {
                coefficients.get(hm_dir.getName()).put(groups_dir.getName(), new HashMap<>());

                File f_data = new File(
                        groups_dir + File.separator + TFPRIO.configs.dynamite.fileStructure.s_output_toBePlotted.get());

                String data = FileManagement.readFile(f_data);

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
        File tf_file = TFPRIO.configs.distributionAnalysis.fileStructure.f_dcg_stats.get();

        File d_plots = TFPRIO.configs.plots.fileStructure.d_output.get();

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
                        String geneID = null;
                        for (String id : TFPRIO.mapSymbolAndEnsg.symbolToEnsg(tf_name))
                        {
                            geneID = id;
                            break;
                        }

                        Map<String, Map<String, Number>> log2fc = new HashMap<>();
                        Map<String, Number> tpm = new HashMap<>();
                        Map<String, Number> normex = new HashMap<>();

                        {   //LOG2FC
                            File d_log2fs = TFPRIO.configs.deSeq2.fileStructure.d_output.get();

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
                            File d_tpm = TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_tpm_tpmResults.get();

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
                            File d_normex = TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_meanCounts.get();


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

        loadTargetGenes();
    }

    private void loadTargetGenes() throws FileNotFoundException
    {
        File tfsDirectory = TFPRIO.configs.distributionAnalysis.fileStructure.d_heatmaps.get();

        for (TranscriptionFactorGroup tfGroup : transcriptionFactorGroups)
        {
            File tfDirectory = FileManagement.getFileIfInDirectory(tfsDirectory, tfGroup.getName(), false);

            if (tfDirectory == null)
            {
                continue;
            }

            for (File d_hm : Objects.requireNonNull(tfDirectory.listFiles()))
            {
                for (File f_data : Objects.requireNonNull(d_hm.listFiles()))
                {
                    if (f_data.getName().endsWith(".csv"))
                    {
                        try (Scanner scanner = new Scanner(f_data))
                        {
                            boolean first = true;
                            int c_geneID = -1, c_geneSymbol = -1;

                            while (scanner.hasNextLine())
                            {
                                String[] line = scanner.nextLine().split(",");

                                if (first)
                                {
                                    for (int i = 0; i < line.length; i++)
                                    {
                                        String strippedEntry = line[i].substring(1, line[i].length() - 1);

                                        if (strippedEntry.equals("geneID"))
                                        {
                                            c_geneID = i;
                                        } else if (strippedEntry.equals("geneSymbol"))
                                        {
                                            c_geneSymbol = i;
                                        }
                                    }

                                    if (c_geneID == -1 || c_geneSymbol == -1)
                                    {
                                        break;
                                    }

                                    first = false;
                                } else
                                {
                                    tfGroup.addTargetGene(line[c_geneSymbol], line[c_geneID]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    public void generate() throws IOException, InterruptedException
    {
        logger.info("Start generating report");

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

        logger.info("Finished generating report");
    }

    private void copyDependencies() throws IOException
    {
        logger.info("Start copying dependencies");
        FileManagement.copyFile(TFPRIO.configs.report.inputStructure.f_style.get(),
                TFPRIO.configs.report.outputStructure.f_style.get());
        FileManagement.copyFile(TFPRIO.configs.report.inputStructure.f_script.get(),
                TFPRIO.configs.report.outputStructure.f_script.get());

        FileManagement.copyDirectory(TFPRIO.configs.report.inputStructure.d_media.get(),
                TFPRIO.configs.report.outputStructure.d_media.get(), false);
        logger.info("Finished copying dependencies");
    }

    private void findExistingValues()
    {
        {
            File pairings = TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_combined.get();

            for (File pairing : Objects.requireNonNull(pairings.listFiles()))
            {
                existingValues.get(SelectorTypes.GROUP_PAIRINGS).add(pairing.getName());
            }
        }

        existingValues.get(SelectorTypes.PERFORMANCE_CUTOFFS).addAll(List.of("1", "2", "3"));

        existingValues.get(SelectorTypes.IMPORTANT_LOCI).addAll(TFPRIO.configs.igv.importantLociAllPrioTf.get());

        existingValues.get(SelectorTypes.TOP_LOG2FC).addAll(Arrays.asList("downregulated", "upregulated"));

        for (Double cutoff : TFPRIO.configs.plots.thresholds.get())
        {
            existingValues.get(SelectorTypes.REGRESSION_CUTOFFS).add(String.valueOf(cutoff));
        }

        existingValues.get(SelectorTypes.DISTRIBUTION_OPTIONS)
                .addAll(existingValues.get(SelectorTypes.HISTONE_MODIFICATIONS));

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
