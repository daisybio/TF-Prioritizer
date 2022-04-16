package lib.AngularReport;

import lib.ExecutableStep;
import org.json.JSONObject;
import tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static util.FileManagement.writeFile;

public class Generate extends ExecutableStep
{
    AbstractConfig<File> f_input_transcriptionFactors = TFPRIO.configs.distributionAnalysis.fileStructure.f_dcg_stats;
    GeneratedFileStructure f_output_reportJson = TFPRIO.configs.angularReport.fileStructure.f_data;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(List.of(f_input_transcriptionFactors));
    }

    @Override public Set<GeneratedFileStructure> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(f_output_reportJson));
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

                TranscriptionFactorGroup transcriptionFactorGroup = new TranscriptionFactorGroup(symbol, logger);
                transcriptionFactorGroups.add(transcriptionFactorGroup);
                executorService.submit(transcriptionFactorGroup::collectData);
            }
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }
        finishAllQueuedThreads();
        saveJson(transcriptionFactorGroups);
    }

    private void saveJson(List<TranscriptionFactorGroup> transcriptionFactorGroups)
    {
        JSONObject jsonObject = new JSONObject();

        for (TranscriptionFactorGroup transcriptionFactorGroup : transcriptionFactorGroups)
        {
            jsonObject.put(transcriptionFactorGroup.name, transcriptionFactorGroup.toJSONObject());
        }

        try
        {
            writeFile(f_output_reportJson.get(), jsonObject.toString(4));
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }
    }
}
