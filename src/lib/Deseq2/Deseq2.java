package lib.Deseq2;

import lib.ExecutableStep;
import tfprio.tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.FileFilters.Filters;

import java.io.*;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;

import static util.FileManagement.*;
import static util.ScriptExecution.executeAndWait;

public class Deseq2 extends ExecutableStep
{
    private final AbstractConfig<File> f_input_scriptDeseq = TFPRIO.configs.scriptTemplates.f_deseq2;
    private final GeneratedFileStructure d_input_combined =
            TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_combined;
    private final GeneratedFileStructure d_output_scripts = TFPRIO.configs.deSeq2.fileStructure.d_rScripts;
    private final GeneratedFileStructure d_output = TFPRIO.configs.deSeq2.fileStructure.d_outputRaw;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(){{
            add(f_input_scriptDeseq);
            add(d_input_combined);
        }};
    }

    @Override public Set<GeneratedFileStructure> getCreatedFileStructure()
    {
        return new HashSet<>(){{
            add(d_output_scripts);
            add(d_output);
        }};
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>();
    }

    @Override protected void execute()
    {

        String inputScriptTemplate = null;
        try
        {
            inputScriptTemplate = readFile(f_input_scriptDeseq.get());
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }
        assert inputScriptTemplate != null;
        String scriptTemplate = inputScriptTemplate;

        for (File f_combination : Objects.requireNonNull(d_input_combined.get().listFiles(Filters.fileFilter)))
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
                    File targetFile = extend(d_output.get(), combination + ".tsv");
                    makeSureFileExists(targetFile);
                    script = script.replace("{OUTPUTFILE}", targetFile.getAbsolutePath());

                    File f_script = extend(d_output_scripts.get(), combination + ".R");
                    writeFile(f_script, script);
                    executeAndWait(f_script, logger);
                } catch (IOException e)
                {
                    logger.error(e.getMessage());
                }
            });
        }
    }
}

