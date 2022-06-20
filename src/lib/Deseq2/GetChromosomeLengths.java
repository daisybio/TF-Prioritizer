package lib.Deseq2;

import lib.ExecutableStep;
import org.json.JSONArray;
import org.json.JSONObject;
import tfprio.tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.IntStream;

import static util.FileManagement.writeFile;

public class GetChromosomeLengths extends ExecutableStep
{
    private final AbstractConfig<String> genomeName = TFPRIO.configs.deSeq2.biomartDatasetSpecies;

    private final GeneratedFileStructure f_output_chromosomeLengths =
            TFPRIO.configs.deSeq2.fileStructure.f_preprocessing_chromosomeLengths;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>();
    }

    @Override public Set<GeneratedFileStructure> getCreatedFileStructure()
    {
        return new HashSet<>()
        {{
            add(f_output_chromosomeLengths);
        }};
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>()
        {{
            add(genomeName);
        }};
    }

    @Override protected void execute()
    {
        StringBuilder response = new StringBuilder();
        try
        {
            URL url = new URL("https://rest.ensembl.org/info/assembly/" + genomeName.get().split("_")[0] + "?");
            HttpURLConnection connection = (HttpURLConnection) url.openConnection();
            connection.setRequestMethod("GET");
            connection.setRequestProperty("Content-Type", "application/json");
            try (BufferedReader reader = new BufferedReader(new InputStreamReader(connection.getInputStream())))
            {
                String inputLine;
                while ((inputLine = reader.readLine()) != null)
                {
                    response.append(inputLine);
                }
            }
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }
        JSONObject json = new JSONObject(response.toString());
        JSONArray topLevelRegions = json.getJSONArray("top_level_region");
        Map<String, Integer> chromosome_length = new HashMap<>();
        IntStream.range(0, topLevelRegions.length()).mapToObj(topLevelRegions::getJSONObject)
                .filter(jsonObject -> jsonObject.getString("coord_system").equals("chromosome"))
                .forEach(entry -> chromosome_length.put(entry.getString("name"), entry.getInt("length")));
        writeFile(f_output_chromosomeLengths.get(), new JSONObject(chromosome_length).toString(4), logger);
    }
}
