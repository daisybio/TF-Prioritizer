package org.exbio.tfprio.lib;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.net.MalformedURLException;
import java.util.*;
import java.util.stream.Collectors;

import static org.exbio.tfprio.lib.Biomart.query;

public class FetchGeneInfo {
    // TODO: Remove retching of sequences and GC content, not used
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

    public FetchGeneInfo(String species, List<String> geneIds) {
        this.geneIds = geneIds;


        List<String[]> exonResponse;
        List<String[]> genesResponse;

        try {
            exonResponse = query(species, geneIds, exonCols);
            genesResponse = query(species, geneIds, geneCols);
        } catch (MalformedURLException e) {
            throw new RuntimeException(e);
        }

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


    public final List<String> getLines() {
        return geneIds.stream().filter(
                id -> geneLengths.containsKey(id) && genePositions.containsKey(id) && geneGc.containsKey(id) &&
                        geneChromosome.containsKey(id) && geneStrand.containsKey(id)).map(
                id -> String.join("\t", id, String.valueOf(geneLengths.get(id)), String.valueOf(geneGc.get(id)),
                        String.valueOf(geneChromosome.get(id)), String.valueOf(genePositions.get(id).getLeft()),
                        String.valueOf(genePositions.get(id).getRight()), String.valueOf(geneStrand.get(id)))).collect(
                Collectors.toList());
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
