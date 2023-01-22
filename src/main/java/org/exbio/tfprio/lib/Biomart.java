package org.exbio.tfprio.lib;

import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLEncoder;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Biomart {
    public static List<String[]> query(String datasetName, List<String> geneIDs, List<String> attributes)
            throws MalformedURLException {
        return query(datasetName, geneIDs, attributes, "ensembl_gene_id");
    }

    public static List<String[]> query(String datasetName, List<String> geneIDs, List<String> attributes,
                                       String geneIdType) throws MalformedURLException {
        List<String> mirrors = List.of("www", "uswest", "useast", "asia");
        String urlString = "http://{MIRROR}.ensembl.org/biomart/martservice";

        String xml = getBiomartXml(datasetName, String.join(",", geneIDs), attributes, geneIdType);

        Map<String, String> parameters = new HashMap<>() {{
            put("query", xml);
        }};

        List<String[]> content = new ArrayList<>();
        boolean successful = false;

        int failedAttempts = 0;

        while (!successful) {
            for (String mirror : mirrors) {
                try {
                    URL url = new URL(urlString.replace("{MIRROR}", mirror));

                    HttpURLConnection connection = (HttpURLConnection) url.openConnection();
                    connection.setRequestMethod("PUT");
                    connection.setDoOutput(true);
                    try (DataOutputStream out = new DataOutputStream(connection.getOutputStream())) {
                        out.writeBytes(getParamsString(parameters));
                        out.flush();
                    }
                    try (BufferedReader reader = new BufferedReader(
                            new InputStreamReader(connection.getInputStream()))) {
                        successful = true;
                        String inputLine;
                        while ((inputLine = reader.readLine()) != null) {
                            content.add(inputLine.split("\t"));
                        }
                    }
                } catch (IOException e) {
                    System.out.println(
                            "Failed to connect to " + mirror + ".ensembl.org. Message was: " + e.getMessage());
                }

                if (successful) {
                    break;
                }
            }

            if (!successful) {
                failedAttempts++;
                if (failedAttempts > 3) {
                    throw new RuntimeException("Failed to connect to biomart");
                }
            }
        }
        return content;
    }

    private static String getBiomartXml(String dataSetName, String geneIDs, List<String> attributes,
                                        String geneIdType) {
        StringBuilder builder = new StringBuilder();
        builder.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        builder.append("<!DOCTYPE Query>\n");
        builder.append(
                "<Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" >\n");
        builder.append("<Dataset name = \"").append(dataSetName).append("\" interface = \"default\" >\n");
        builder.append("<Filter name = \"").append(geneIdType).append("\" value = \"").append(geneIDs).append(
                "\" />\n");
        builder.append("<Attribute name = \"").append(geneIdType).append("\" />\n");
        for (String attribute : attributes) {
            builder.append("<Attribute name = \"").append(attribute).append("\" />\n");
        }
        builder.append("</Dataset>\n");
        builder.append("</Query>\n");
        return builder.toString();
    }

    private static String getParamsString(Map<String, String> params) {
        StringBuilder builder = new StringBuilder();

        for (Map.Entry<String, String> entry : params.entrySet()) {
            builder.append(URLEncoder.encode(entry.getKey(), StandardCharsets.UTF_8));
            builder.append("=");
            builder.append(URLEncoder.encode(entry.getValue(), StandardCharsets.UTF_8));
            builder.append("&");
        }

        String string = builder.toString();
        return string.length() > 0 ? string.substring(0, string.length() - 1) : string;
    }
}
