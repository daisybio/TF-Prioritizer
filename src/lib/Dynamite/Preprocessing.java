package lib.Dynamite;

import lib.ExecutableStep;
import tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.FileManagement;

import java.io.File;
import java.util.*;

import static util.FileManagement.extend;
import static util.FileManagement.makeSureFileExists;
import static util.ScriptExecution.executeAndWait;

public class Preprocessing extends ExecutableStep
{
    private AbstractConfig<File> d_input;
    private final AbstractConfig<File> d_inputDeseq2 = TFPRIO.configs.deSeq2.fileStructure.d_output;
    private final AbstractConfig<File> f_script = TFPRIO.configs.tepic.f_dynamite_integrateDate;
    private final AbstractConfig<File> f_scriptClassification =
            TFPRIO.configs.tepic.f_dynamite_prepareForClassification;

    private final GeneratedFileStructure d_output = TFPRIO.configs.dynamite.fileStructure.d_preprocessing_integrateData;
    private final GeneratedFileStructure d_outputClassification =
            TFPRIO.configs.dynamite.fileStructure.d_preprocessing_prepareClassification;

    private final AbstractConfig<String> s_ratios_dir =
            TFPRIO.configs.tepic.fileStructure.s_postprocessing_output_ratiosDir;
    private final AbstractConfig<String> s_log2coeff =
            TFPRIO.configs.dynamite.fileStructure.s_preprocessing_integrateData_log2coeff;
    private final AbstractConfig<String> s_prepClass =
            TFPRIO.configs.dynamite.fileStructure.s_preprocessing_prepareClassification_data;
    private final AbstractConfig<Integer> geneIdColumn = TFPRIO.configs.dynamite.preprocessing_IntegrateDataGeneIds;
    private final AbstractConfig<Integer> log2fcColumn = TFPRIO.configs.dynamite.preprocessing_IntegrateDataLog2fc;

    private final AbstractConfig<File> tgeneExecutable = TFPRIO.configs.tgene.pathToExecutable;
    private final AbstractConfig<File> considerGene =
            TFPRIO.configs.dynamite.preprocessing_IntegrateDataConsiderGeneFile;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(Arrays.asList(d_input, d_inputDeseq2, f_script, f_scriptClassification));
    }

    @Override public Set<GeneratedFileStructure> getCreatedFileStructure()
    {
        return new HashSet<>(Arrays.asList(d_output, d_outputClassification));
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>(Arrays.asList(s_ratios_dir, s_log2coeff, s_prepClass, geneIdColumn, log2fcColumn));
    }

    @Override protected Set<AbstractConfig<?>> getOptionalConfigs()
    {
        return new HashSet<>(Arrays.asList(tgeneExecutable, considerGene));
    }

    @Override protected void updateInputDirectory()
    {
        d_input = FileManagement.getFirstExisting(Arrays.asList(TFPRIO.configs.tgene.fileStructure.d_integrate,
                TFPRIO.configs.tgene.fileStructure.d_filteredTargetGenes,
                TFPRIO.configs.tepic.fileStructure.d_postprocessing_output));
        logger.debug("Input: " + d_input.get().getAbsolutePath());
    }

    @Override protected void execute()
    {
        for (Map.Entry<String, Set<String>> pairingEntry : TFPRIO.groupCombinationsToHms.entrySet())
        {
            String pairing = pairingEntry.getKey();

            for (String hm : pairingEntry.getValue())
            {
                File input_ratios;

                if (!tgeneExecutable.isSet())
                {
                    input_ratios = extend(d_input.get(), pairing, hm, s_ratios_dir.get(), pairing + ".txt");
                } else
                {
                    input_ratios = extend(d_input.get(), pairing, hm, pairing + ".txt");
                }

                File input_diff_gene_expr = extend(d_inputDeseq2.get(), pairing + ".tsv");
                File targetFile = extend(d_output.get(), pairing, hm, s_log2coeff.get());

                makeSureFileExists(targetFile, logger);

                List<String> parameters =
                        Arrays.asList(input_ratios.getAbsolutePath(), input_diff_gene_expr.getAbsolutePath(),
                                targetFile.getAbsolutePath(), "--geneIDs " + geneIdColumn.get(),
                                "--expressionC " + log2fcColumn.get());
                if (considerGene.isSet())
                {
                    parameters.add("--filterIDs " + considerGene.get().getAbsolutePath());
                }

                String command = "python3 " + f_script.get().getAbsolutePath() + " " + String.join(" ", parameters);

                executorService.submit(() -> executeAndWait(command, logger));
            }
        }

        finishAllQueuedThreads();

        for (Map.Entry<String, Set<String>> pairingEntry : TFPRIO.groupCombinationsToHms.entrySet())
        {
            String pairing = pairingEntry.getKey();

            for (String hm : pairingEntry.getValue())
            {
                File targetFile = extend(d_outputClassification.get(), pairing, hm, s_prepClass.get());
                makeSureFileExists(targetFile, logger);

                String command = "Rscript " + f_scriptClassification.get().getAbsolutePath() + " " +
                        extend(d_output.get(), pairing, hm, s_log2coeff.get()).getAbsolutePath() + " " +
                        targetFile.getAbsolutePath();

                executorService.submit(() -> executeAndWait(command, logger));
            }
        }
    }
}
