package org.exbio.tfprio.lib;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.net.MalformedURLException;
import java.util.*;
import java.util.stream.Collectors;

import static org.exbio.tfprio.lib.Biomart.query;

public class FetchGeneInfo {
    private static final List<String> exonCols = List.of("ensembl_exon_id", "exon_chrom_start", "exon_chrom_end");
    private static final List<String> geneCols = List.of("chromosome_name", "start_position", "end_position");
    private static final Logger logger = LogManager.getLogger("Biomart");
    private final List<String> geneIds;
    private final Map<String, Integer> geneLengths = new HashMap<>();
    private final Map<String, Pair<Integer, Integer>> genePositions = new HashMap<>();
    private final Map<String, String> geneChromosome = new HashMap<>();

    public FetchGeneInfo(String species, List<String> geneIds, String geneIdType) {
        this.geneIds = geneIds;

        List<String[]> exonResponse;
        List<String[]> genesResponse;

        try {
            exonResponse = query(species, geneIds, exonCols, geneIdType);
            genesResponse = query(species, geneIds, geneCols, geneIdType);
        } catch (MalformedURLException e) {
            throw new RuntimeException(e);
        }

        Map<String, List<String[]>> exonMap = mapIds(exonResponse);
        Map<String, List<String[]>> genesMap = mapIds(genesResponse);

        geneIds.forEach(id -> {
            List<String[]> exonLines = exonMap.getOrDefault(id, new ArrayList<>());
            List<String[]> geneLines = genesMap.getOrDefault(id, new ArrayList<>());

            if (geneLines.isEmpty()) {
                return;
            }

            if (geneLines.size() > 1) {
                logger.warn("Gene {} has {} gc lines! Only processing the first one.", id, geneLines.size());
            }

            String[] geneLine = geneLines.get(0);

            if (geneLine.length != (geneCols.size() + 1)) {
                logger.warn("Gene {} has {} columns!", id, geneLine.length);
                logger.debug(Arrays.toString(geneLine));
            }

            if (exonLines.isEmpty()) {
                return;
            }

            int cumulatedExonLengths = exonLines.stream().mapToInt(line -> {
                int start = Integer.parseInt(line[exonCols.indexOf("exon_chrom_start") + 1]);
                int end = Integer.parseInt(line[exonCols.indexOf("exon_chrom_end") + 1]);
                return Math.abs(end - start);
            }).sum();

            geneLengths.put(id, cumulatedExonLengths);

            int start = Integer.parseInt(geneLine[geneCols.indexOf("start_position") + 1]);
            int end = Integer.parseInt(geneLine[geneCols.indexOf("end_position") + 1]);
            String chromosome = geneLine[geneCols.indexOf("chromosome_name") + 1];

            genePositions.put(id, Pair.of(start, end));
            geneChromosome.put(id, chromosome);
        });
    }


    public final List<String> getLines() {
        return geneIds.stream().filter(id -> geneLengths.containsKey(id) && genePositions.containsKey(id) &&
                geneChromosome.containsKey(id)).map(
                id -> String.join("\t", id, String.valueOf(geneLengths.get(id)), String.valueOf(geneChromosome.get(id)),
                        String.valueOf(genePositions.get(id).getLeft()),
                        String.valueOf(genePositions.get(id).getRight()))).collect(Collectors.toList());
    }

    private Map<String, List<String[]>> mapIds(List<String[]> matrix) {
        return new HashMap<>() {{
            for (String[] entry : matrix) {
                String id = entry[0];
                computeIfAbsent(id, $ -> new ArrayList<>()).add(entry);
            }
        }};
    }
}
