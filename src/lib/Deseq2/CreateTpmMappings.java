package lib.Deseq2;

import lib.ExecutableStep;
import tfprio.TFPRIO;
import util.Configs.Config;
import util.FileFilters.Filters;
import util.Hashing;

import java.io.*;
import java.net.HttpURLConnection;
import java.net.URL;
import java.net.URLEncoder;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.regex.Pattern;

import static util.FileManagement.*;
import static util.ScriptExecution.executeAndWait;

public class CreateTpmMappings extends ExecutableStep
{
    private final Config<File> d_meanCounts = TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_meanCounts;
    private final Config<File> d_output = TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_tpm_tpmResults;
    private final Config<File> d_outputScripts = TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_tpm_scripts;
    private final Config<File> f_lengths = TFPRIO.configs.deSeq2.fileStructure.f_preprocessing_tpm_geneLengths;
    private final Config<File> f_geneIDs = TFPRIO.configs.deSeq2.inputGeneID;

    private final Config<File> f_scriptTemplate = TFPRIO.configs.scriptTemplates.f_deseq2PreprocessingTpm;

    private final Config<Boolean> calculateTpmLengthsEnabled = TFPRIO.configs.general.calculateTpmLengthsEnabled;
    private final Config<Integer> threadLimit = TFPRIO.configs.general.threadLimit;
    private final Config<String> species = TFPRIO.configs.deSeq2.biomartDatasetSpecies;

    @Override protected Set<Config<File>> getRequiredFileStructure()
    {
        return new HashSet<>(Arrays.asList(d_meanCounts, f_geneIDs, f_scriptTemplate));
    }

    @Override protected Set<Config<File>> getCreatedFileStructure()
    {
        return new HashSet<>(Arrays.asList(d_output, d_outputScripts, f_lengths));
    }

    @Override protected Set<Config<?>> getRequiredConfigs()
    {
        return new HashSet<>(Arrays.asList(calculateTpmLengthsEnabled, threadLimit, species));
    }


    private String getBiomartXml(String dataSetName, String geneIDs, List<String> attributes)
    {
        StringBuilder body = new StringBuilder("<?xml version='1.0' encoding='UTF-8'?><!DOCTYPE Query>");
        body.append("<Query  virtualSchemaName = 'default' formatter='TSV' header='0' uniqueRows='0' count='' " +
                "datasetConfigVersion='0.6'>");
        body.append("<Dataset name='").append(dataSetName).append("' interface='default'>");
        body.append("<Filter name='ensembl_gene_id' value ='").append(geneIDs).append("'/>");

        for (String attribute : attributes)
        {
            body.append("<Attribute name='").append(attribute).append("'/>");
        }

        body.append("</Dataset>");
        body.append("</Query>");

        return body.toString();
    }

    public static String getParamsString(Map<String, String> params)
    {
        StringBuilder result = new StringBuilder();

        for (Map.Entry<String, String> entry : params.entrySet())
        {
            result.append(URLEncoder.encode(entry.getKey(), StandardCharsets.UTF_8));
            result.append("=");
            result.append(URLEncoder.encode(entry.getValue(), StandardCharsets.UTF_8));
            result.append("&");
        }

        String resultString = result.toString();
        return resultString.length() > 0 ? resultString.substring(0, resultString.length() - 1) : resultString;
    }

    private List<String[]> query(String datasetName, List<String> geneIDs, List<String> attributes)
            throws IOException, InterruptedException
    {
        String urlString = "http://www.ensembl.org/biomart/martservice";

        StringBuilder sb_geneID = new StringBuilder();
        for (String geneID : geneIDs)
        {
            sb_geneID.append(geneID).append(",");
        }
        sb_geneID.deleteCharAt(sb_geneID.lastIndexOf(","));
        String xml = getBiomartXml(datasetName, sb_geneID.toString(), attributes);
        URL url = new URL(urlString);


        Map<String, String> parameters = new HashMap<>()
        {{
            put("query", xml);
        }};

        List<String[]> content = new ArrayList<>();
        boolean successful = false;

        while (!successful)
        {
            HttpURLConnection connection = (HttpURLConnection) url.openConnection();
            connection.setRequestMethod("PUT");
            connection.setDoOutput(true);
            try (DataOutputStream out = new DataOutputStream(connection.getOutputStream()))
            {
                out.writeBytes(getParamsString(parameters));
                out.flush();
            }
            try (BufferedReader reader = new BufferedReader(new InputStreamReader(connection.getInputStream())))
            {
                successful = true;
                String inputLine;
                while ((inputLine = reader.readLine()) != null)
                {
                    content.add(inputLine.split("\t"));
                }
            } catch (IOException e)
            {
                logger.warn("Retrying: " + Hashing.hash(String.valueOf(geneIDs.hashCode())));
            }
        }
        return content;
    }

    private List<String[]> reduce(List<String[]> input, List<String> columns)
    {
        Set<Integer> addedIndices = new HashSet<>();
        List<String[]> output = new ArrayList<>();

        for (int i = 0; i < input.size(); i++)
        {
            if (addedIndices.contains(i))
            {
                continue;
            }
            addedIndices.add(i);

            String[] entry = input.get(i);
            int start = Integer.parseInt(entry[columns.indexOf("exon_chrom_start")]);
            int end = Integer.parseInt(entry[columns.indexOf("exon_chrom_end")]);

            for (int j = 0; j < input.size(); j++)
            {
                if (addedIndices.contains(j))
                {
                    continue;
                }
                String[] compareEntry = input.get(j);
                int compareStart = Integer.parseInt(compareEntry[columns.indexOf("exon_chrom_start")]);
                int compareEnd = Integer.parseInt(compareEntry[columns.indexOf("exon_chrom_end")]);

                if (compareEnd < start || compareStart > end)
                {
                    continue;
                }

                start = Math.min(start, compareStart);
                end = Math.max(end, compareEnd);
                addedIndices.add(j);
            }

            entry[columns.indexOf("exon_chrom_start")] = String.valueOf(start);
            entry[columns.indexOf("exon_chrom_end")] = String.valueOf(end);
            output.add(entry);
        }

        return output;
    }

    private Map<String, List<String[]>> mapIds(List<String> columns, List<String[]> matrix)
    {
        Map<String, List<String[]>> map = new HashMap<>();
        for (String[] entry : matrix)
        {
            String id = entry[columns.indexOf("ensembl_gene_id")];
            if (!map.containsKey(id))
            {
                map.put(id, new ArrayList<>());
            }
            map.get(id).add(entry);
        }
        return map;
    }

    @Override protected void execute()
    {
        logger.info("Create TPM values for all RNA-seq data.");

        logger.info("Get gene lengths ...");

        List<String> lengthCols =
                Arrays.asList("ensembl_gene_id", "ensembl_exon_id", "chromosome_name", "exon_chrom_start",
                        "exon_chrom_end");
        List<String> gcCols = Arrays.asList("gene_exon_intron", "ensembl_gene_id", "start_position", "end_position");
        try
        {
            List<String> geneIDs = readLines(f_geneIDs.get());
            geneIDs.remove(0);
            Map<String, Map<String, Number>> result = new HashMap<>();

            int splitNumber = threadLimit.get() * 2;
            int splitSize = (geneIDs.size() + splitNumber - 1) / splitNumber;
            logger.info("GeneIDs: " + geneIDs.size() + ", batches: " + splitNumber + ", batchSize: " + splitSize);

            for (int splitIndex = 0; splitIndex < splitNumber; splitIndex++)
            {
                int finalSplitIndex = splitIndex;
                executorService.submit(() ->
                {
                    int startIndex = finalSplitIndex * splitSize;
                    int endIndex = Math.min((finalSplitIndex + 1) * splitSize, geneIDs.size());
                    List<String> selectedIDs = geneIDs.subList(startIndex, endIndex);
                    List<String[]> r_length = null;
                    List<String[]> r_gc = null;
                    try
                    {
                        r_length = query(species.get(), selectedIDs, lengthCols);
                        r_gc = query(species.get(), selectedIDs, gcCols);
                    } catch (IOException | InterruptedException e)
                    {
                        e.printStackTrace();
                    }
                    assert r_length != null;
                    assert r_gc != null;

                    Map<String, List<String[]>> m_length = mapIds(lengthCols, r_length);
                    Map<String, List<String[]>> m_gc = mapIds(gcCols, r_gc);

                    for (String id : m_length.keySet())
                    {
                        List<String[]> l_length = m_length.get(id);
                        List<String[]> l_length_reduced = reduce(l_length, lengthCols);
                        result.put(id, new HashMap<>());

                        int length = 0;
                        for (String[] entry : l_length_reduced)
                        {
                            length += Integer.parseInt(entry[lengthCols.indexOf("exon_chrom_end")]) -
                                    Integer.parseInt(entry[lengthCols.indexOf("exon_chrom_start")]) + 1;
                        }
                        result.get(id).put("Length", length);

                        if (m_gc.containsKey(id))
                        {
                            List<String[]> l_gc = m_gc.get(id);

                            assert l_gc.size() == 1;

                            String sequence = l_gc.get(0)[gcCols.indexOf("gene_exon_intron")];

                            int offset = Integer.MAX_VALUE;
                            for (String[] entry : l_length_reduced)
                            {
                                int start = Integer.parseInt(entry[lengthCols.indexOf("exon_chrom_start")]);
                                offset = Math.min(start, offset);
                            }
                            StringBuilder sb_exonSequence = new StringBuilder();
                            for (String[] entry : l_length_reduced)
                            {
                                int start = Integer.parseInt(entry[lengthCols.indexOf("exon_chrom_start")]);
                                int end = Integer.parseInt(entry[lengthCols.indexOf("exon_chrom_end")]);

                                sb_exonSequence.append(sequence, start - offset, end - offset + 1);
                            }

                            String exonSequence = sb_exonSequence.toString();

                            double gcContent =
                                    (double) Pattern.compile("[GC]").matcher(exonSequence).results().count() /
                                            exonSequence.length();

                            result.get(id).put("gcContent", gcContent);
                        }
                    }
                    logger.debug("Batch finished: " + startIndex + "-" + endIndex);
                });
            }

            finishAllQueuedThreads();

            makeSureFileExists(f_lengths.get());
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(f_lengths.get())))
            {
                writer.write("ENSG\tlength\tgc");
                writer.newLine();
                for (String geneID : geneIDs)
                {
                    writer.write(geneID);
                    writer.write("\t");
                    if (result.containsKey(geneID))
                    {
                        writer.write(String.valueOf(result.get(geneID).get("Length")));
                        writer.write("\t");
                        if (result.get(geneID).containsKey("gcContent"))
                        {
                            writer.write(String.valueOf(result.get(geneID).get("gcContent")));
                        } else
                        {
                            writer.write("NA");
                        }
                    } else
                    {
                        writer.write("NA");
                        writer.write("\t");
                        writer.write("NA");
                    }
                    writer.newLine();
                }
            }

        } catch (IOException e)
        {
            e.printStackTrace();
        }


        try
        {
            String scriptTemplate = readFile(f_scriptTemplate.get());

            for (File f_group : Objects.requireNonNull(d_meanCounts.get().listFiles(Filters.fileFilter)))
            {
                executorService.submit(() ->
                {
                    String group = f_group.getName().substring(0, f_group.getName().lastIndexOf("."));
                    String script = scriptTemplate.replace("{COUNTFILE}", f_group.getAbsolutePath());
                    script = script.replace("{LENGTHSFILE}", f_lengths.get().getAbsolutePath());
                    File targetFile = extend(d_output.get(), group + ".tsv");
                    try
                    {
                        makeSureFileExists(targetFile);
                        script = script.replace("{TARGETFILE}", targetFile.getAbsolutePath());

                        File scriptFile = extend(d_outputScripts.get(), group + ".py");

                        writeFile(scriptFile, script);
                        executeAndWait(scriptFile, logger);
                    } catch (IOException e)
                    {
                        logger.error(e.getMessage());
                        System.exit(1);
                    }
                });
            }
        } catch (IOException e)
        {
            e.printStackTrace();
            System.exit(1);
        }
    }
}
