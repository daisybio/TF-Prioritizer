package lib;

import tfprio.TFPRIO;
import util.ExternalScriptException;
import util.FileFilters.Filters;

import java.io.*;
import java.util.*;

import static util.FileManagement.*;
import static util.ScriptExecution.executeAndWait;

public class CreateDeseq2Scripts extends ExecutableStep
{
    @Override protected void execute()
    {
        File d_outputCombined = TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_combined.get();
        File d_outputSingle = TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_single.get();
        File d_outputMeanCounts = TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_meanCounts.get();

        File input = TFPRIO.configs.deSeq2.inputDirectory.get();


        // Single
        try
        {
            String scriptTemplate = readFile(TFPRIO.configs.scriptTemplates.f_deseq2PreprocessingSingle.get());
            for (File d_group : Objects.requireNonNull(input.listFiles(Filters.directoryFilter)))
            {
                File d_outputGroup = extend(d_outputSingle, d_group.getName());
                String script = scriptTemplate.replace("{INPUT_DIRECTORY}", d_group.getAbsolutePath());
                script = script.replace("{OUTPUT_DIRECTORY}", d_outputGroup.getAbsolutePath());
                File meansFile = extend(d_outputMeanCounts, d_group.getName() + ".tsv");
                makeSureFileExists(meansFile);
                script = script.replace("{OUTPUT_MEANS_FILE}", meansFile.getAbsolutePath());
                script = script.replace("{SYMBOL_MAP_FILE}",
                        TFPRIO.configs.deSeq2.fileStructure.f_mapping.get().getAbsolutePath());
                String finalScript =
                        script.replace("{GENE_ID_FILE}", TFPRIO.configs.deSeq2.inputGeneID.get().getAbsolutePath());

                executorService.execute(() ->
                {
                    try
                    {
                        executeAndWait(finalScript, ".py");
                    } catch (ExternalScriptException e)
                    {
                        logger.error(e.getMessage());
                        System.exit(1);
                    }
                });
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
            String scriptTemplate = readFile(TFPRIO.configs.scriptTemplates.f_deseq2PreprocessingCombined.get());
            for (File d_group1 : Objects.requireNonNull(d_outputSingle.listFiles(Filters.directoryFilter)))
            {
                for (File d_group2 : Objects.requireNonNull(d_outputSingle.listFiles(Filters.directoryFilter)))
                {
                    if (d_group1.getName().compareTo(d_group2.getName()) >= 0)
                    {
                        continue;
                    }
                    String combination = d_group1.getName() + "_" + d_group2.getName();
                    String script = scriptTemplate.replace("{INPUT_DIRECTORY1}", d_group1.getAbsolutePath());
                    script = script.replace("{INPUT_DIRECTORY2}", d_group2.getAbsolutePath());
                    script = script.replace("{GENE_ID}", TFPRIO.configs.deSeq2.inputGeneID.get().getAbsolutePath());
                    File targetFile = extend(d_outputCombined, combination + ".tsv");
                    makeSureFileExists(targetFile);
                    String finalScript = script.replace("{OUTPUT_FILE}", targetFile.getAbsolutePath());

                    executorService.execute(() ->
                    {
                        try
                        {
                            executeAndWait(finalScript, ".py");
                        } catch (ExternalScriptException e)
                        {
                            logger.error(e.getMessage());
                            System.exit(1);
                        }
                    });
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
        File r_scripts = TFPRIO.configs.deSeq2.fileStructure.d_rScripts.get();

        try
        {
            String scriptTemplate = readFile(TFPRIO.configs.scriptTemplates.f_deseq2.get());
            for (File f_combination : Objects.requireNonNull(d_outputCombined.listFiles(Filters.fileFilter)))
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
                        for (String sample : samples)
                        {
                            String group = sample.substring(0, sample.indexOf("_"));
                            sb_samples.append("'" + sample + "', ");
                            sb_groups.append("'" + group + "', ");
                        }
                        sb_samples.setLength(sb_samples.length() - 2);
                        sb_groups.setLength(sb_groups.length() - 2);
                        script = script.replace("{ SAMPLES }", sb_samples);
                        script = script.replace("{ GROUPS }", sb_groups);
                        script = script.replace("{INPUTFILE}", f_combination.getAbsolutePath());
                        File targetFile =
                                extend(TFPRIO.configs.deSeq2.fileStructure.d_outputRaw.get(), combination + ".tsv");
                        makeSureFileExists(targetFile);
                        script = script.replace("{OUTPUTFILE}", targetFile.getAbsolutePath());

                        writeFile(extend(r_scripts, combination + ".R"), script);
                    } catch (IOException e)
                    {
                        e.printStackTrace();
                        System.exit(1);
                    }
                });
            }
        } catch (IOException e)
        {
            e.printStackTrace();
        }

    }
}
