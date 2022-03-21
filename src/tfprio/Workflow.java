package tfprio;

import lib.ExecutableStep;
import util.Logger;

import java.util.ArrayList;
import java.util.List;

public class Workflow
{
    private final List<ExecutableStep> steps = new ArrayList<>();
    private final Logger logger = new Logger("Workflow");

    public Workflow()
    {
        steps.add(new tfprio.InitStaticVariables());
        steps.add(new lib.CheckChromosomes());

        if (TFPRIO.configs.mixOptions.mutuallyExclusive.get() && !TFPRIO.configs.mixOptions.level.isSet())
        {
            try
            {
                TFPRIO.configs.mixOptions.level.setValue("SAMPLE_LEVEL");
            } catch (IllegalAccessException ignore)
            {
            }
        }

        if (TFPRIO.configs.mixOptions.level.isSet())
        {
            steps.add(new lib.MixOptions.SampleLevelPreprocess());
            steps.add(new lib.MixOptions.SampleLevel());

            if (TFPRIO.configs.mixOptions.level.get().equals("HM_LEVEL"))
            {
                steps.add(new lib.MixOptions.HmLevelPreprocess());
                steps.add(new lib.MixOptions.HmLevel());
            }
        }

        if (TFPRIO.configs.tepic.tfBindingSiteSearch.isSet() &&
                (TFPRIO.configs.tepic.tfBindingSiteSearch.get().equals("BETWEEN") ||
                        TFPRIO.configs.tepic.tfBindingSiteSearch.get().equals("EXCL_BETWEEN")))
        {
            steps.add(new lib.FootprintIntervals.CreateFootprintsBetweenPeaks());
        }

        if (TFPRIO.configs.blacklist.bedFilePath.isSet())
        {
            steps.add(new lib.Blacklist.Blacklist());
        }

        if (TFPRIO.configs.mixOptions.mutuallyExclusive.get())
        {
            steps.add(new lib.MixOptions.MixMutuallyExclusive());
        }

        {
            steps.add(new lib.Deseq2.CreateDeseq2Scripts());
            steps.add(new lib.Deseq2.CreateTpmMappings());
            steps.add(new lib.Deseq2.CreateGenePositions());

            if (TFPRIO.configs.deSeq2.tpmFilter.isSet())
            {
                steps.add(new lib.Deseq2.PreprocessTpm());
            }
            steps.add(new lib.Deseq2.Deseq2());
            steps.add(new lib.Deseq2.PostProcessing());
        }

        if (TFPRIO.configs.tgene.pathToExecutable.isSet())
        {
            steps.add(new lib.Tgene.RunTgene());
            steps.add(new lib.Tgene.Postprocessing());
        }
    }

    public boolean simulationSuccessful()
    {
        boolean allWorked = true;
        for (ExecutableStep step : steps)
        {
            allWorked = allWorked && step.simulate();
        }
        return allWorked;
    }

    public void filterByHashes()
    {
        long startTime = System.currentTimeMillis();
        logger.info("Start filtering by hashes.");

        ArrayList<ExecutableStep> filtered = new ArrayList<>();
        boolean previousFailed = false;

        for (ExecutableStep step : steps)
        {
            if (!previousFailed)
            {
                if (step.verifyHash())
                {
                    logger.info("Found valid hash for: " + step.getClass().getName());
                } else
                {
                    logger.info("Hash for " + step.getClass().getName() +
                            " is invalid. Not using hashed files from now on.");
                    filtered.add(step);
                    previousFailed = true;
                }
            } else
            {
                filtered.add(step);
            }
        }

        steps.clear();
        steps.addAll(filtered);

        double deltaSeconds = (double) (System.currentTimeMillis() - startTime) / 1e3;
        logger.info("Finished filtering by hashes. Execution took " + deltaSeconds + " seconds.");
    }

    public void createHashes()
    {
        long startTime = System.currentTimeMillis();
        logger.info("Start creating hashes.");
        for (ExecutableStep step : steps)
        {
            step.createHash();
        }
        double deltaSeconds = (double) (System.currentTimeMillis() - startTime) / 1e3;
        logger.info("Finished creating hashes. Execution took " + deltaSeconds + " seconds.");
    }

    public void run()
    {
        if (TFPRIO.configs.general.hashingEnabled.get())
        {
            filterByHashes();
        }

        logger.info("Start running executable steps.");
        for (ExecutableStep step : steps)
        {
            step.run();
        }
        logger.info("Finished running executable steps.");

        if (TFPRIO.configs.general.hashingEnabled.get())
        {
            createHashes();
        }
    }

    public String toString()
    {
        StringBuilder sb_out = new StringBuilder();

        for (ExecutableStep step : steps)
        {
            sb_out.append(step.getClass().getName()).append("\n");
        }
        return sb_out.toString();
    }
}
