package lib.Igv;

import lib.ExecutableStep;
import tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.FileFilters.Filters;
import util.IGV_Headless;

import java.io.File;
import java.util.*;

import static lib.Igv.Helpers.getGeneCoordinates;
import static lib.Igv.Helpers.getInputFiles;
import static tfprio.Workflow.getLatestInputDirectory;
import static util.FileManagement.extend;

public class ImportantLoci extends ExecutableStep
{
    private final AbstractConfig<File> f_input_geneCoordinates =
            TFPRIO.configs.deSeq2.fileStructure.f_preprocessing_genePositions_data;
    private AbstractConfig<File> d_input_tepic;
    private final AbstractConfig<File> d_input_peakFiles = TFPRIO.configs.chipAtlas.fileStructure.d_peakFiles;
    private final AbstractConfig<File> d_input_dcg =
            TFPRIO.configs.distributionAnalysis.fileStructure.d_dcg_targetGenes;

    private final GeneratedFileStructure d_output = TFPRIO.configs.igv.fileStructure.d_importantLoci;

    private final AbstractConfig<List<String>> importantLoci = TFPRIO.configs.igv.importantLociAllPrioTf;
    private final AbstractConfig<String> s_session = TFPRIO.configs.igv.fileStructure.s_session;
    private final AbstractConfig<String> speciesReferenceGenome = TFPRIO.configs.igv.speciesReferenceGenome;

    private final AbstractConfig<List<String>> includePredictionData = TFPRIO.configs.igv.includePredictionData;
    private final AbstractConfig<File> pathToTfChipSeq = TFPRIO.configs.igv.pathToTfChipSeq;
    private final AbstractConfig<File> pathToTdf = TFPRIO.configs.igv.pathToTdf;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(Arrays.asList(f_input_geneCoordinates, d_input_tepic, d_input_peakFiles, d_input_dcg));
    }

    @Override public Set<GeneratedFileStructure> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(d_output));
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>(Arrays.asList(importantLoci, s_session, speciesReferenceGenome));
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
        Map<String, String> gene_coordinates = getGeneCoordinates(f_input_geneCoordinates.get(), logger);

        ArrayList<String> add_to_important_loci = new ArrayList<>();
        ArrayList<String> remove_from_important_loci = new ArrayList<>();

        //check if symbol is available, if not check if similar one is there
        for (Object locusObject : importantLoci.get())
        {
            String locus = (String) locusObject;
            if (!gene_coordinates.containsKey(locus))
            {
                for (String loci : gene_coordinates.keySet())
                {
                    if (loci.contains(locus))
                    {
                        add_to_important_loci.add(loci);
                    }
                }
                remove_from_important_loci.add(locus);
            }
        }
        importantLoci.get().addAll(add_to_important_loci);
        importantLoci.get().removeAll(remove_from_important_loci);

        for (File d_group : Objects.requireNonNull(d_input_dcg.get().listFiles(Filters.directoryFilter)))
        {
            executorService.submit(() ->
            {
                String group = d_group.getName();

                File d_output_group = extend(d_output.get(), group);

                List<File> tdfFiles = new ArrayList<>();

                List<String> loadFiles =
                        getInputFiles(List.of(group), includePredictionData, d_input_tepic, pathToTfChipSeq, pathToTdf,
                                d_input_peakFiles, tdfFiles);

                String loadCommand = "load " + String.join("\nload ", loadFiles);

                File f_save_session = extend(d_output_group, s_session.get());
                IGV_Headless.createSession(f_save_session, loadCommand, tdfFiles, logger);

                IGV_Headless igv = new IGV_Headless(group, logger);

                igv.addCommand("genome " + speciesReferenceGenome.get());
                igv.addCommand("load " + f_save_session.getAbsolutePath());
                igv.addCommand("snapshotDirectory " + d_output_group.getAbsolutePath());

                for (Object locusObject : importantLoci.get())
                {
                    String locus = (String) locusObject;

                    if (!gene_coordinates.containsKey(locus))
                    {
                        continue;
                    }

                    igv.addCommand("goto " + gene_coordinates.get(locus));
                    igv.addCommand("snapshot " + locus + ".png");
                }

                igv.run(new File(d_output_group.getAbsolutePath()));
            });
        }
    }
}
