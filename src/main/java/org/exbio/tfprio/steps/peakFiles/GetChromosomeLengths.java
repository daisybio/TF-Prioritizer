package org.exbio.tfprio.steps.peakFiles;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;
import org.json.JSONArray;
import org.json.JSONObject;

import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.IntStream;

public class GetChromosomeLengths extends ExecutableStep {
    public final OutputFile outputFile = addOutput("chromosomeLengths.json");
    private final RequiredConfig<String> species = new RequiredConfig<>(Configs.mixOptions.biomartDatasetSpecies);

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                StringBuilder response = new StringBuilder();

                // Fetch data from the ensembl API
                try {
                    URL url = new URL("https://rest.ensembl.org/info/assembly/" + species.get().split("_")[0] + "?");
                    HttpURLConnection connection = (HttpURLConnection) url.openConnection();
                    connection.setRequestMethod("GET");
                    connection.setRequestProperty("Content-Type", "application/json");

                    try (BufferedReader reader = new BufferedReader(
                            new InputStreamReader(connection.getInputStream()))) {
                        String inputLine;
                        while ((inputLine = reader.readLine()) != null) {
                            response.append(inputLine);
                        }
                    }
                } catch (IOException e) {
                    logger.error(e.getMessage());
                }

                // Parse the data
                JSONObject json = new JSONObject(response.toString());
                JSONArray topLevelRegions = json.getJSONArray("top_level_region");
                Map<String, Integer> chromosome_length = new HashMap<>();
                IntStream.range(0, topLevelRegions.length()).mapToObj(topLevelRegions::getJSONObject).filter(
                        jsonObject -> jsonObject.getString("coord_system").equals("chromosome")).forEach(
                        entry -> chromosome_length.put(entry.getString("name"), entry.getInt("length")));

                // Write the data to a file
                try (FileWriter writer = new FileWriter(outputFile)) {
                    writer.write(new JSONObject(chromosome_length).toString(4));
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
                return true;
            });
        }};
    }
}
