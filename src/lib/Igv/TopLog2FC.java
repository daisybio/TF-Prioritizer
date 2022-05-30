package lib.Igv;

import lib.ExecutableStep;
import lib.GeneAffinity;
import tfprio.tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.FileFilters.Filters;
import util.IGV_Headless;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

import static lib.Igv.Helpers.*;
import static tfprio.tfprio.Workflow.getLatestInputDirectory;
import static util.FileManagement.extend;

public class TopLog2FC extends ExecutableStep
{
    private final AbstractConfig<File> f_input_geneCoordinates =
            TFPRIO.configs.deSeq2.fileStructure.f_preprocessing_genePositions_data;
    private final AbstractConfig<File> d_input_deseq2 = TFPRIO.configs.deSeq2.fileStructure.d_output;
    private AbstractConfig<File> d_input_tepic;
    private final AbstractConfig<File> d_input_peakFiles = TFPRIO.configs.chipAtlas.fileStructure.d_peakFiles;
    private final AbstractConfig<File> d_input_bedFiles =
            TFPRIO.configs.tepic.fileStructure.d_postprocessing_trapPredictedBeds;

    private final GeneratedFileStructure d_output = TFPRIO.configs.igv.fileStructure.d_igvTopLog2fc;

    private final AbstractConfig<String> speciesReferenceGenome = TFPRIO.configs.igv.speciesReferenceGenome;
    private final AbstractConfig<String> s_session = TFPRIO.configs.igv.fileStructure.s_session;
    private final AbstractConfig<String> s_upregulated = TFPRIO.configs.igv.fileStructure.s_igvTopLog2fc_upregulated;
    private final AbstractConfig<String> s_downregulated =
            TFPRIO.configs.igv.fileStructure.s_igvTopLog2fc_downregulated;
    private final AbstractConfig<Integer> topKLog2FCs = TFPRIO.configs.igv.topLog2fc;
    private final AbstractConfig<Boolean> includeLncRnaPseudogenes =
            TFPRIO.configs.igv.topLog2fcIncludeLncRnaPseudogenes;

    private final AbstractConfig<List<String>> includePredictionData = TFPRIO.configs.igv.includePredictionData;
    private final AbstractConfig<File> pathToTfChipSeq = TFPRIO.configs.igv.pathToTfChipSeq;
    private final AbstractConfig<File> pathToTdf = TFPRIO.configs.igv.pathToTdf;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(Arrays.asList(f_input_geneCoordinates, d_input_deseq2, d_input_tepic, d_input_peakFiles,
                d_input_peakFiles));
    }

    @Override public Set<GeneratedFileStructure> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(d_output));
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>(
                Arrays.asList(speciesReferenceGenome, s_session, s_upregulated, s_downregulated, topKLog2FCs,
                        includeLncRnaPseudogenes));
    }

    @Override protected Set<AbstractConfig<?>> getOptionalConfigs()
    {
        return new HashSet<>(Arrays.asList(includePredictionData, pathToTfChipSeq, pathToTdf));
    }

    @Override protected void updateInputDirectory()
    {
        d_input_tepic = getLatestInputDirectory();
    }

    @Override protected void execute()
    {
        List<String> geneClasses = Arrays.asList("lncRNA", "pseudogene", "miRNA");
        Map<String, String> gene_coordinates = getGeneCoordinates(f_input_geneCoordinates.get(), logger);

        HashMap<String, String> ensg_class = new HashMap<>();

        for (String gene : gene_coordinates.keySet())
        {
            try
            {
                String description = TFPRIO.mapSymbolAndEnsg.anyToDescription(gene);
                for (String geneClass : geneClasses)
                {
                    if (description.contains(geneClass))
                    {
                        Set<String> ensgs = new HashSet<>();

                        if (TFPRIO.mapSymbolAndEnsg.hasSymbol(gene))
                        {
                            ensgs.addAll(TFPRIO.mapSymbolAndEnsg.symbolToEnsg(gene));
                        } else if (TFPRIO.mapSymbolAndEnsg.hasEnsg(gene))
                        {
                            ensgs.add(gene);
                        }
                        for (String ensg : ensgs)
                        {
                            ensg_class.put(ensg, geneClass);
                        }
                        break;
                    }
                }
            } catch (NoSuchFieldException e)
            {
                logger.warn(e.getMessage());
            }
        }

        for (File f_groupPairing : Objects.requireNonNull(d_input_deseq2.get().listFiles(Filters.fileFilter)))
        {
            executorService.submit(() ->
            {
                String groupPairing = f_groupPairing.getName().substring(0, f_groupPairing.getName().lastIndexOf("."));

                ArrayList<GeneAffinity> geneAffinities = new ArrayList<>();
                try (BufferedReader br_input_log2fc = new BufferedReader(new FileReader(f_groupPairing)))
                {
                    String line_input_log2fc;
                    br_input_log2fc.readLine();
                    while ((line_input_log2fc = br_input_log2fc.readLine()) != null)
                    {
                        String[] split = line_input_log2fc.split("\t");
                        GeneAffinity geneAffinity = new GeneAffinity(split[0], Double.parseDouble(split[1]));
                        geneAffinities.add(geneAffinity);
                    }
                } catch (IOException e)
                {
                    logger.error(e.getMessage());
                }
                Collections.sort(geneAffinities);

                File f_save_session = extend(d_output.get(), groupPairing, s_session.get());
                ArrayList<File> tdfFiles = new ArrayList<>();

                List<String> loadFiles =
                        getInputFiles(List.of(groupPairing.split("_")), includePredictionData, d_input_tepic,
                                pathToTfChipSeq, pathToTdf, d_input_peakFiles, tdfFiles);
                addBedFiles(loadFiles, List.of(groupPairing.split("_")), TFPRIO.existingHms, d_input_bedFiles);

                for (String suffix : Arrays.asList(s_downregulated.get(), s_upregulated.get()))
                {
                    File d_output_mode = extend(d_output.get(), groupPairing, suffix);
                    ArrayList<GeneAffinity> regulatedGenes = new ArrayList<>();

                    for (int i = 0; i < topKLog2FCs.get(); i++)
                    {
                        final GeneAffinity geneAffinity = geneAffinities.get(i);
                        final String geneClass;

                        geneClass = ensg_class.getOrDefault(geneAffinity.getGeneID(), "");

                        if (TFPRIO.mapSymbolAndEnsg.hasEnsg(geneAffinity.getGeneID()) || includeLncRnaPseudogenes.get())
                        {
                            if (!geneClass.isEmpty())
                            {
                                geneAffinity.setGeneClass(geneClass);
                            }
                            regulatedGenes.add(geneAffinity);
                        }
                    }

                    IGV_Headless igv = new IGV_Headless(groupPairing + "-" + suffix, logger);
                    igv.createSession(loadFiles, tdfFiles, f_save_session);

                    igv.addCommand("snapshotDirectory " + d_output_mode.getAbsolutePath());

                    for (int i = 0; i < regulatedGenes.size(); i++)
                    {
                        int rank = i + 1;
                        String locus = regulatedGenes.get(i).getGeneSymbol();

                        if (!gene_coordinates.containsKey(locus))
                        {
                            continue;
                        }

                        String snapshot_name = rank + "_" + regulatedGenes.get(i).getGeneSymbol() + ".png";

                        igv.addCommand("goto " + gene_coordinates.get(locus));
                        igv.addCommand("snapshot " + snapshot_name);
                    }

                    igv.run(d_output_mode);
                }
            });
        }
    }
}
