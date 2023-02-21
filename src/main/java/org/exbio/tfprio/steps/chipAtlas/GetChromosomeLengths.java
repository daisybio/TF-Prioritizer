package org.exbio.tfprio.steps.chipAtlas;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;
import org.json.JSONArray;
import org.json.JSONObject;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.Collection;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class GetChromosomeLengths extends ExecutableStep<Configs> {
    public final OutputFile outputFile = addOutput("chromosome_lengths.tsv");
    private final RequiredConfig<String> species = new RequiredConfig<>(configs.inputConfigs.biomartSpecies);

    public GetChromosomeLengths(Configs configs) {
        super(configs);
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                URL url = new URL("https://rest.ensembl.org/info/assembly/" + species.get().split("_")[0] + "?");
                HttpURLConnection connection = (HttpURLConnection) url.openConnection();
                connection.setRequestMethod("GET");
                connection.setRequestProperty("Content-Type", "application/json");

                final JSONObject response;

                try (BufferedReader reader = new BufferedReader(new InputStreamReader(connection.getInputStream()))) {
                    response = new JSONObject(reader.lines().collect(Collectors.joining()));
                }

                JSONArray topLevelRegions = response.getJSONArray("top_level_region");

                Map<String, Integer> chromosomeLengths =
                        IntStream.range(0, topLevelRegions.length()).mapToObj(topLevelRegions::getJSONObject).filter(
                                obj -> obj.getString("coord_system").equals("chromosome")).collect(
                                Collectors.toMap(entry -> entry.getString("name"), entry -> entry.getInt("length")));

                try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
                    for (Map.Entry<String, Integer> entry : chromosomeLengths.entrySet()) {
                        writer.write(entry.getKey() + "\t" + entry.getValue());
                        writer.newLine();
                    }
                }

                return true;
            });
        }};
    }
}
