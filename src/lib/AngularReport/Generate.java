package lib.AngularReport;

import lib.ExecutableStep;
import org.json.JSONObject;
import tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.Logger;
import util.ScriptExecution;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ExecutorService;

import static util.FileManagement.*;
import static util.ScriptExecution.executeAndWait;

public class Generate extends ExecutableStep
{
    AbstractConfig<File> f_input_transcriptionFactors = TFPRIO.configs.distributionAnalysis.fileStructure.f_dcg_stats;
    GeneratedFileStructure f_output_reportJson = TFPRIO.configs.angularReport.fileStructure.f_data;
    GeneratedFileStructure d_output_data = TFPRIO.configs.angularReport.fileStructure.d_data;


    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(List.of(f_input_transcriptionFactors));
    }

    @Override public Set<GeneratedFileStructure> getCreatedFileStructure()
    {
        return new HashSet<>(Arrays.asList(f_output_reportJson, d_output_data));
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>();
    }

    @Override protected void execute()
    {
        List<TranscriptionFactorGroup> transcriptionFactorGroups = new ArrayList<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(f_input_transcriptionFactors.get())))
        {
            String inputLine;
            reader.readLine();

            while ((inputLine = reader.readLine()) != null)
            {
                String[] split = inputLine.split("\t");

                if (split.length != 3)
                {
                    logger.warn("Irregular line length in tf ranking file detected");
                }

                String symbol = split[1];

                TranscriptionFactorGroup transcriptionFactorGroup =
                        new TranscriptionFactorGroup(symbol, logger, executorService);
                transcriptionFactorGroups.add(transcriptionFactorGroup);
                transcriptionFactorGroup.collectData();
            }
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }
        finishAllQueuedThreads();
        saveJson(transcriptionFactorGroups);

        softLink(TFPRIO.configs.angularReport.d_dataLink.get(),
                TFPRIO.configs.angularReport.fileStructure.d_preprocessing.get(), logger);

        String script =
                "cd " + TFPRIO.configs.angularReport.d_angularSource.get().getAbsolutePath() + "\n" + "ng build " +
                        "--output-path=" + TFPRIO.configs.angularReport.fileStructure.d_output.get().getAbsolutePath();

        try
        {
            writeFile(TFPRIO.configs.angularReport.fileStructure.f_script.get(), script);
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        executeAndWait(TFPRIO.configs.angularReport.fileStructure.f_script.get(), logger);
    }

    private void saveJson(List<TranscriptionFactorGroup> transcriptionFactorGroups)
    {
        JSONObject jsonObject = new JSONObject();

        List<JSONObject> tfGroupObjects = new ArrayList<>();

        for (TranscriptionFactorGroup transcriptionFactorGroup : transcriptionFactorGroups)
        {
            tfGroupObjects.add(transcriptionFactorGroup.toJSONObject());
        }

        jsonObject.put("configs", TFPRIO.configs.getConfigsJSONObject(true));
        jsonObject.put("transcriptionFactorGroups", tfGroupObjects);

        try
        {
            writeFile(f_output_reportJson.get(), jsonObject.toString(4));
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }
    }

    public static void linkFiles(JSONObject map, File d_target, ExecutorService executorService, Logger logger)
            throws IOException
    {
        for (String key : map.keySet())
        {
            File newTarget = extend(d_target, key);

            Object value = map.get(key);

            if (value.getClass().equals(String.class))
            {
                File sourceFile = new File((String) value);
                String extension = sourceFile.getName().substring(sourceFile.getName().lastIndexOf("."));

                File targetFile = new File(newTarget.getAbsolutePath() + extension);

                executorService.submit(() -> softLink(targetFile, sourceFile, logger));

                map.put(key, targetFile.getAbsolutePath()
                        .replace(TFPRIO.configs.angularReport.fileStructure.d_data.get().getAbsolutePath(), ""));

            } else if (value.getClass().equals(JSONObject.class))
            {
                linkFiles((JSONObject) value, newTarget, executorService, logger);
            } else
            {
                throw new IllegalArgumentException("Illegal map structure detected");
            }
        }
    }
}