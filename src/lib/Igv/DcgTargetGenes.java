package lib.Igv;

import lib.ExecutableStep;
import tfprio.TFPRIO;
import tfprio.Workflow;
import util.Configs.Config;
import util.FileFilters.Filters;
import util.IGV_Headless;

import java.io.*;
import java.util.*;

import static util.FileManagement.extend;

public class DcgTargetGenes extends ExecutableStep
{
    private final Config<File> f_input_geneCoordinates =
            TFPRIO.configs.deSeq2.fileStructure.f_preprocessing_genePositions_data;
    private final Config<File> d_input_heatmaps = TFPRIO.configs.distributionAnalysis.fileStructure.d_heatmaps;
    private Config<File> d_input_tepic;
    private final Config<File> d_input_peakFiles = TFPRIO.configs.chipAtlas.fileStructure.d_peakFiles;

    private final Config<File> d_output = TFPRIO.configs.igv.fileStructure.d_igvDcgTargetGenes;

    private final Config<String> speciesReferenceGenome = TFPRIO.configs.igv.speciesReferenceGenome;
    private final Config<String> s_session = TFPRIO.configs.igv.fileStructure.s_session;
    private final Config<File> pathToIgv = TFPRIO.configs.igv.pathToIGV;
    private final Config<Boolean> chipAtlasEnabled = TFPRIO.configs.chipAtlas.isEnabled;

    // Optional configs
    private final Config<List> includePredictionData = TFPRIO.configs.igv.includePredictionData;
    private final Config<File> pathToTfChipSeq = TFPRIO.configs.igv.pathToTfChipSeq;
    private final Config<File> pathToTdf = TFPRIO.configs.igv.pathToTdf;
    private final Config<List> enhancerDatabases = TFPRIO.configs.igv.enhancerDatabases;

    @Override protected Set<Config<File>> getRequiredFileStructure()
    {
        return new HashSet<>(Arrays.asList(f_input_geneCoordinates, d_input_heatmaps, d_input_tepic))
        {{
            if (TFPRIO.configs.igv.enhancerDatabases.isSet() && TFPRIO.configs.igv.enhancerDatabases.get().size() > 0)
            {
                add(TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_genePositions_enhancerDBs);
            }

            if (chipAtlasEnabled.get())
            {
                add(d_input_peakFiles);
            }
        }};
    }

    @Override protected Set<Config<File>> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(d_output));
    }

    @Override protected Set<Config<?>> getRequiredConfigs()
    {
        return new HashSet<>(Arrays.asList(speciesReferenceGenome, s_session, pathToIgv, chipAtlasEnabled));
    }

    @Override protected Set<Config<?>> getOptionalConfigs()
    {
        return new HashSet<>(Arrays.asList(includePredictionData, pathToTfChipSeq, pathToTdf, enhancerDatabases));
    }

    @Override protected void updateInputDirectory()
    {
        d_input_tepic = Workflow.getLatestInputDirectory();
    }

    @Override protected void execute()
    {
        HashMap<String, String> gene_coordinates = new HashMap<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(f_input_geneCoordinates.get())))
        {
            String inputLine;
            reader.readLine();

            while ((inputLine = reader.readLine()) != null)
            {
                String[] split = inputLine.split("\t");

                if (split[1].startsWith("CHR"))
                {
                    continue;
                }

                String gene_name;
                if (split[0].equals(""))
                {
                    gene_name = split[6];
                } else
                {
                    gene_name = split[0];
                }

                String chr = "chr" + split[1];
                int begin = Integer.parseInt(split[2]);
                int end = Integer.parseInt(split[3]);

                begin -= 50000;
                end += 50000;

                String position_string = chr + ":" + begin + "-" + end;
                gene_coordinates.put(gene_name.toUpperCase(), position_string);
            }
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }


        for (File d_tf : Objects.requireNonNull(d_input_heatmaps.get().listFiles(Filters.directoryFilter)))
        {
            String tfSymbol = d_tf.getName();

            for (File d_hm : Objects.requireNonNull(d_tf.listFiles(Filters.directoryFilter)))
            {
                String hm = d_hm.getName();

                for (File f_groupPairing : Objects.requireNonNull(d_hm.listFiles(Filters.getSuffixFilter(".csv"))))
                {
                    executorService.submit(() ->
                    {


                        String groupPairing =
                                f_groupPairing.getName().substring(0, f_groupPairing.getName().lastIndexOf("."));

                        String[] groupsSplit = groupPairing.split("_");

                        ArrayList<String> regulated_genes = new ArrayList<>();
                        //get gene names to get shots of
                        try (BufferedReader br_genes = new BufferedReader(new FileReader(f_groupPairing)))
                        {
                            String line_genes;
                            br_genes.readLine();
                            while ((line_genes = br_genes.readLine()) != null)
                            {
                                String[] split = line_genes.split(",");
                                String ensg = split[0].replace("\"", "");

                                if (TFPRIO.mapSymbolAndEnsg.hasEnsg(ensg))
                                {
                                    try
                                    {
                                        regulated_genes.add(TFPRIO.mapSymbolAndEnsg.ensgToSymbol(ensg));
                                    } catch (NoSuchFieldException ignore)
                                    {
                                    }
                                }
                            }
                        } catch (IOException e)
                        {
                            logger.error(e.getMessage());
                        }

                        File d_output_groupPairing = extend(d_output.get(), tfSymbol, hm, groupPairing);

                        //include all peak file etc in IGV
                        //create IGV load for both groups
                        ArrayList<File> tdf_files = new ArrayList<>();

                        ArrayList<String> loadFiles = new ArrayList<>();

                        for (String group : groupsSplit)
                        {
                            //include predictive HMs
                            if (includePredictionData.isSet() && includePredictionData.get().size() > 0)
                            {
                                for (Object includingHmObject : includePredictionData.get())
                                {
                                    String includingHm = (String) includingHmObject;

                                    File d_input = extend(d_input_tepic.get(), group, includingHm);

                                    if (d_input.exists())
                                    {
                                        for (File f_input : Objects.requireNonNull(
                                                d_input.listFiles(Filters.fileFilter)))
                                        {
                                            loadFiles.add(f_input.getAbsolutePath());
                                        }
                                    }
                                }
                            }

                            //add all own TF ChIP-seq if available
                            if (pathToTfChipSeq.isSet())
                            {
                                File d_input = extend(pathToTfChipSeq.get(), group);
                                if (d_input.exists() && d_input.isDirectory())
                                {
                                    for (File d_input_tf : Objects.requireNonNull(
                                            d_input.listFiles(Filters.directoryFilter)))
                                    {
                                        for (File f_input : Objects.requireNonNull(
                                                d_input_tf.listFiles(Filters.fileFilter)))
                                        {
                                            loadFiles.add(f_input.getAbsolutePath());
                                        }
                                    }
                                }
                            }


                            //add TDF if available
                            if (pathToTdf.isSet())
                            {
                                File d_input = extend(pathToTdf.get(), group);
                                if (d_input.exists() && d_input.isDirectory())
                                {
                                    for (File d_input_hm : Objects.requireNonNull(
                                            d_input.listFiles(Filters.directoryFilter)))
                                    {
                                        for (File f_input : Objects.requireNonNull(
                                                d_input_hm.listFiles(Filters.fileFilter)))
                                        {
                                            loadFiles.add(f_input.getAbsolutePath());
                                            tdf_files.add(f_input);
                                        }
                                    }
                                }
                            }
                        }

                        if (d_input_peakFiles.get().exists())
                        {
                            for (File d_input_tf : Objects.requireNonNull(
                                    d_input_peakFiles.get().listFiles(Filters.directoryFilter)))
                            {
                                for (File f_input : Objects.requireNonNull(d_input_tf.listFiles(Filters.fileFilter)))
                                {
                                    String file_name = f_input.getName();

                                    if (file_name.endsWith("idx") || file_name.endsWith("bam"))
                                    {
                                        continue;
                                    }

                                    loadFiles.add(f_input.getAbsolutePath());
                                }
                            }
                        }

                        String loadCommand = "load " + String.join(", ", loadFiles);

                        //save session
                        File f_save_session = extend(d_output_groupPairing, s_session.get());

                        IGV_Headless.createSession(f_save_session, loadCommand, tdf_files, logger);

                        IGV_Headless igv = new IGV_Headless(groupPairing, logger);

                        igv.addCommand("genome" + speciesReferenceGenome);
                        igv.addCommand("load " + f_save_session.getAbsolutePath());

                        for (int i = 0; i < regulated_genes.size(); i++)
                        {
                            int rank = i + 1;
                            String locus = regulated_genes.get(i);

                            if (!gene_coordinates.containsKey(locus))
                            {
                                continue;
                            }

                            igv.addCommand("snapshotDirectory " + d_output_groupPairing.getAbsolutePath());
                            igv.addCommand("goto " + gene_coordinates.get(locus));
                            igv.addCommand("snapshot " + rank + "_" + regulated_genes.get(i) + ".png");
                        }
                        igv.run(d_output_groupPairing);
                    });
                }
            }
        }
    }
}
