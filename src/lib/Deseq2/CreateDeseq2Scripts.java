package lib.Deseq2;

import lib.ExecutableStep;
import tfprio.tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.FileFilters.Filters;

import java.io.*;
import java.util.*;

import static util.FileManagement.*;
import static util.ScriptExecution.executeAndWait;

public class CreateDeseq2Scripts extends ExecutableStep
{
    private final AbstractConfig<File> d_input = TFPRIO.configs.deSeq2.inputDirectory;
    private final AbstractConfig<File> f_batches = TFPRIO.configs.deSeq2.batchFile;
    private final GeneratedFileStructure d_outputCombined =
            TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_combined;
    private final GeneratedFileStructure d_outputSingle = TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_single;
    private final GeneratedFileStructure d_outputMeanCounts =
            TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_meanCounts;
    private final GeneratedFileStructure d_outputScripts = TFPRIO.configs.deSeq2.fileStructure.d_rScripts;
    private final GeneratedFileStructure d_outputRaw = TFPRIO.configs.deSeq2.fileStructure.d_outputRaw;

    private final AbstractConfig<File> f_scriptSingle = TFPRIO.configs.scriptTemplates.f_deseq2PreprocessingSingle;
    private final AbstractConfig<File> f_scriptCombined = TFPRIO.configs.scriptTemplates.f_deseq2PreprocessingCombined;
    private final AbstractConfig<File> f_scriptDeseq = TFPRIO.configs.scriptTemplates.f_deseq2;

    private final AbstractConfig<File> f_geneID = TFPRIO.configs.deSeq2.inputGeneID;
    private final AbstractConfig<File> f_mapping = TFPRIO.configs.deSeq2.fileStructure.f_mapping;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(
                Arrays.asList(d_input, f_scriptCombined, f_scriptDeseq, f_scriptSingle, f_geneID, f_mapping))
        {{
            if (f_batches.isSet())
            {
                add(f_batches);
            }
        }};
    }

    @Override public Set<GeneratedFileStructure> getCreatedFileStructure()
    {
        return new HashSet<>(
                Arrays.asList(d_outputSingle, d_outputCombined, d_outputMeanCounts, d_outputScripts, d_outputRaw));
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>();
    }

    @Override protected void execute()
    {
        // Single
        try
        {
            String scriptTemplate = readFile(f_scriptSingle.get());

            for (File d_group : Objects.requireNonNull(d_input.get().listFiles(Filters.directoryFilter)))
            {
                File d_outputGroup = extend(d_outputSingle.get(), d_group.getName());
                String script = scriptTemplate.replace("{INPUT_DIRECTORY}", d_group.getAbsolutePath());
                script = script.replace("{OUTPUT_DIRECTORY}", d_outputGroup.getAbsolutePath());
                File meansFile = extend(d_outputMeanCounts.get(), d_group.getName() + ".tsv");
                makeSureFileExists(meansFile);
                script = script.replace("{OUTPUT_MEANS_FILE}", meansFile.getAbsolutePath());
                script = script.replace("{SYMBOL_MAP_FILE}", f_mapping.get().getAbsolutePath());
                String finalScript = script.replace("{GENE_ID_FILE}", f_geneID.get().getAbsolutePath());

                executorService.execute(() -> executeAndWait(finalScript, ".py", logger));
            }
        } catch (IOException e)
        {
            e.printStackTrace();
            System.exit(1);
        }

        finishAllQueuedThreads();

        // Combined
        try
        {
            String scriptTemplate = readFile(f_scriptCombined.get());
            for (File d_group1 : Objects.requireNonNull(d_outputSingle.get().listFiles(Filters.directoryFilter)))
            {
                for (File d_group2 : Objects.requireNonNull(d_outputSingle.get().listFiles(Filters.directoryFilter)))
                {
                    if (d_group1.getName().compareTo(d_group2.getName()) >= 0)
                    {
                        continue;
                    }
                    String combination = d_group1.getName() + "_" + d_group2.getName();
                    String script = scriptTemplate.replace("{INPUT_DIRECTORY1}", d_group1.getAbsolutePath());
                    script = script.replace("{INPUT_DIRECTORY2}", d_group2.getAbsolutePath());
                    script = script.replace("{GENE_ID}", f_geneID.get().getAbsolutePath());
                    File targetFile = extend(d_outputCombined.get(), combination + ".tsv");
                    makeSureFileExists(targetFile);
                    String finalScript = script.replace("{OUTPUT_FILE}", targetFile.getAbsolutePath());

                    executorService.execute(() -> executeAndWait(finalScript, ".py", logger));
                }
            }
        } catch (IOException e)
        {
            e.printStackTrace();
            System.exit(1);
        }

        finishAllQueuedThreads();

        logger.info("Start creating RScripts for running DESeq2");

        // Scripts
        try
        {
            String scriptTemplate = readFile(f_scriptDeseq.get());
            for (File f_combination : Objects.requireNonNull(d_outputCombined.get().listFiles(Filters.fileFilter)))
            {
                executorService.execute(() ->
                {
                    String combination = f_combination.getName().substring(0, f_combination.getName().lastIndexOf("."));
                    String script = scriptTemplate.replace("{COMBINATION}", combination);
                    String firstLine;
                    try (BufferedReader reader = new BufferedReader(new FileReader(f_combination)))
                    {
                        firstLine = reader.readLine();
                        String[] samples = firstLine.substring(firstLine.indexOf("\t") + 1).split("\t");

                        StringBuilder sb_samples = new StringBuilder();
                        StringBuilder sb_groups = new StringBuilder();
                        StringBuilder sb_batches = new StringBuilder(", batch = c(");

                        for (String sample : samples)
                        {
                            sb_samples.append("'").append(sample.replace("-", ".")).append("', ");

                            if (TFPRIO.sample_group.containsKey(sample) &&
                                    (TFPRIO.sample_batch == null || TFPRIO.sample_batch.containsKey(sample)))
                            {
                                sb_groups.append("'").append(TFPRIO.sample_group.get(sample)).append("', ");
                                if (TFPRIO.sample_batch != null)
                                {
                                    sb_batches.append("'").append(TFPRIO.sample_batch.get(sample)).append("', ");
                                }
                            } else
                            {
                                logger.error("Unknown sample: " + sample);
                            }
                        }
                        sb_samples.setLength(sb_samples.length() - 2);
                        sb_groups.setLength(sb_groups.length() - 2);
                        if (TFPRIO.sample_batch != null)
                        {
                            sb_batches.setLength(sb_batches.length() - 2);
                            sb_batches.append(")");
                            script = script.replace("{ BATCHES }", sb_batches);
                        } else
                        {
                            script = script.replace("{ BATCHES }", "");
                        }
                        script = script.replace("{ DESIGN }", TFPRIO.sample_batch != null ? "batch+group" : "group");
                        script = script.replace("{ SAMPLES }", sb_samples);
                        script = script.replace("{ GROUPS }", sb_groups);
                        script = script.replace("{INPUTFILE}", f_combination.getAbsolutePath());
                        File targetFile = extend(d_outputRaw.get(), combination + ".tsv");
                        makeSureFileExists(targetFile);
                        script = script.replace("{OUTPUTFILE}", targetFile.getAbsolutePath());

                        writeFile(extend(d_outputScripts.get(), combination + ".R"), script);
                    } catch (IOException e)
                    {
                        logger.error(e.getMessage());
                    }
                });
            }
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

    }
}
