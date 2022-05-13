package lib.Dynamite;

import lib.ExecutableStep;
import tfprio.tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.FileFilters.Filters;

import java.io.*;
import java.util.*;

import static util.FileManagement.*;
import static util.ScriptExecution.executeAndWait;

public class RunDynamite extends ExecutableStep
{
    private final AbstractConfig<File> f_dynamite_script = TFPRIO.configs.tepic.f_dynamite;
    private final AbstractConfig<File> d_input =
            TFPRIO.configs.dynamite.fileStructure.d_preprocessing_prepareClassification;

    private final GeneratedFileStructure d_output = TFPRIO.configs.dynamite.fileStructure.d_output;

    private final AbstractConfig<String> outVar = TFPRIO.configs.dynamite.outVar;
    private final AbstractConfig<Integer> oFolds = TFPRIO.configs.dynamite.oFolds;
    private final AbstractConfig<Integer> iFolds = TFPRIO.configs.dynamite.iFolds;
    private final AbstractConfig<Boolean> performance = TFPRIO.configs.dynamite.performance;
    private final AbstractConfig<Double> alpha = TFPRIO.configs.dynamite.alpha;
    private final AbstractConfig<Integer> cores = TFPRIO.configs.dynamite.cores;
    private final AbstractConfig<Boolean> randomize = TFPRIO.configs.dynamite.randomize;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(Arrays.asList(f_dynamite_script, d_input));
    }

    @Override public Set<GeneratedFileStructure> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(d_output));
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>(Arrays.asList(outVar, oFolds, iFolds, performance, alpha, cores, randomize));
    }

    @Override protected void execute()
    {
        for (Map.Entry<String, Set<String>> pairingEntry : TFPRIO.groupCombinationsToHms.entrySet())
        {
            String pairing = pairingEntry.getKey();

            for (String hm : pairingEntry.getValue())
            {
                File d_input_pairingHm = extend(d_input.get(), pairing, hm);

                for (File f_input : Objects.requireNonNull(d_input_pairingHm.listFiles(Filters.fileFilter)))
                {
                    try
                    {
                        if (readLines(f_input).size() <= 2)
                        {
                            logger.error("File corrupted: "  + f_input.getAbsolutePath());
                        }
                    } catch (IOException e)
                    {
                        logger.error(e.getMessage());
                    }
                }

                File d_out = extend(d_output.get(), pairing, hm);
                makeSureDirectoryExists(d_out, logger);

                String command_edited = "Rscript " + f_dynamite_script.get().getAbsolutePath() + " --dataDir=" +
                        d_input_pairingHm.getAbsolutePath() + " --outDir=" + d_out.getAbsolutePath() + " --out_var=" + outVar.get() +
                        " --Ofolds=" + oFolds.get() + " --Ifolds=" + iFolds.get() + " --performance=" +
                        performance.get() + " --alpha=" + alpha.get() + " --cores=" + cores.get() + " --randomise=" +
                        randomize.get();

                executorService.submit(() -> executeAndWait(command_edited, logger));
            }
        }
    }
}
