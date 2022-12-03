package org.exbio.tfprio.lib;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.exbio.pipejar.util.Hashing;

import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLEncoder;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.stream.Collectors;

public class Biomart {
    private static final List<String> exonCols =
            List.of("ensembl_gene_id", "ensembl_exon_id", "chromosome_name", "exon_chrom_start", "exon_chrom_end");
    private static final List<String> geneCols =
            List.of("gene_exon_intron", "ensembl_gene_id", "chromosome_name", "start_position", "end_position",
                    "strand");
    private static final Logger logger = LogManager.getLogger("Biomart");
    private final List<String> geneIds;
    private final Map<String, Integer> geneLengths = new HashMap<>();
    private final Map<String, Double> geneGc = new HashMap<>();
    private final Map<String, Pair<Integer, Integer>> genePositions = new HashMap<>();
    private final Map<String, String> geneChromosome = new HashMap<>();
    private final Map<String, Integer> geneStrand = new HashMap<>();

    public Biomart(String species, List<String> geneIds) {
        this.geneIds = geneIds;

        List<String[]> exonResponse = query(species, geneIds, exonCols);
        List<String[]> genesResponse = query(species, geneIds, geneCols);

        Map<String, List<String[]>> exonMap = mapIds(exonCols, exonResponse);
        Map<String, List<String[]>> genesMap = mapIds(geneCols, genesResponse);

        geneIds.forEach(id -> {
            List<String[]> exonLines = exonMap.getOrDefault(id, new ArrayList<>());
            List<String[]> geneLines = genesMap.getOrDefault(id, new ArrayList<>());

            int offset = Integer.MAX_VALUE;
            Set<Pair<Integer, Integer>> exonCoords = new HashSet<>();

            if (!exonLines.isEmpty()) {
                int length = 0;

                for (String[] exonLine : exonLines) {

                    int start = Integer.parseInt(exonLine[exonCols.indexOf("exon_chrom_start")]);
                    int end = Integer.parseInt(exonLine[exonCols.indexOf("exon_chrom_end")]);

                    int smaller = Math.min(start, end);
                    int larger = Math.max(start, end);

                    length += larger - smaller;

                    exonCoords.add(Pair.of(smaller, larger));

                    offset = Math.min(offset, smaller);
                }

                geneLengths.put(id, length);
            }

            if (geneLines.isEmpty()) {
                return;
            }

            if (geneLines.size() > 1) {
                logger.warn("Gene {} has {} gc lines! Only processing the first one.", id, geneLines.size());
            }

            String[] geneLine = geneLines.get(0);

            if (geneLine.length != geneCols.size()) {
                logger.warn("Gene {} has {} columns!", id, geneLine.length);
                System.out.println(Arrays.toString(geneLine));
            }

            String sequence = geneLine[geneCols.indexOf("gene_exon_intron")];

            int start = Integer.parseInt(geneLine[geneCols.indexOf("start_position")]);
            int end = Integer.parseInt(geneLine[geneCols.indexOf("end_position")]);
            String chromosome = geneLine[geneCols.indexOf("chromosome_name")];
            int strand = Integer.parseInt(geneLine[geneCols.indexOf("strand")]);

            genePositions.put(id, Pair.of(start, end));
            geneChromosome.put(id, chromosome);
            geneStrand.put(id, strand);

            final int finalOffset = offset;

            String exonSequence = exonCoords.stream().map(
                    coord -> sequence.substring(coord.getLeft() - finalOffset, coord.getRight() - finalOffset)).collect(
                    Collectors.joining());
            geneGc.put(id,
                    exonSequence.chars().filter(c -> c == 'G' || c == 'C').count() / (double) exonSequence.length());
        });
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

    private static String getBiomartXml(String dataSetName, String geneIDs, List<String> attributes) {
        StringBuilder builder = new StringBuilder();
        builder.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        builder.append("<!DOCTYPE Query>\n");
        builder.append(
                "<Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" >\n");
        builder.append("<Dataset name = \"").append(dataSetName).append("\" interface = \"default\" >\n");
        builder.append("<Filter name = \"ensembl_gene_id\" value = \"").append(geneIDs).append("\" />\n");
        for (String attribute : attributes) {
            builder.append("<Attribute name = \"").append(attribute).append("\" />\n");
        }
        builder.append("</Dataset>\n");
        builder.append("</Query>\n");
        return builder.toString();
    }

    private static List<String[]> query(String datasetName, List<String> geneIDs, List<String> attributes) {
        String urlString = "http://www.ensembl.org/biomart/martservice";

        String xml = getBiomartXml(datasetName, String.join(",", geneIDs), attributes);

        try {
            URL url = new URL(urlString);
            Map<String, String> parameters = new HashMap<>() {{
                put("query", xml);
            }};

            List<String[]> content = new ArrayList<>();
            boolean successful = false;

            while (!successful) {
                try {
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
                    logger.warn(e.getMessage());
                    logger.warn("Retrying: " + Hashing.hash(String.valueOf(geneIDs.hashCode())));
                }
            }
            return content;
        } catch (MalformedURLException e) {
            logger.error(e.getMessage());
        }
        return new ArrayList<>();
    }

    public final List<String> getLines() {
        return geneIds.stream().map(id -> String.join("\t", id, String.valueOf(geneLengths.getOrDefault(id, null)),
                String.valueOf(geneGc.getOrDefault(id, null)), String.valueOf(geneChromosome.get(id)),
                String.valueOf(genePositions.containsKey(id) ? genePositions.get(id).getLeft() : null),
                String.valueOf(genePositions.containsKey(id) ? genePositions.get(id).getRight() : null),
                String.valueOf(geneStrand.get(id)))).collect(Collectors.toList());
    }

    private Map<String, List<String[]>> mapIds(List<String> columns, List<String[]> matrix) {
        return new HashMap<>() {{
            for (String[] entry : matrix) {
                String id = entry[columns.indexOf("ensembl_gene_id")];
                computeIfAbsent(id, $ -> new ArrayList<>()).add(entry);
            }
        }};
    }
}
