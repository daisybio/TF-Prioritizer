package tfprio;

import lib.ExecutableStep;
import util.ExecutionTimeMeasurement;
import util.Logger;

import java.util.ArrayList;
import java.util.List;

public class Workflow
{
    private final List<ExecutableStep> steps = new ArrayList<>();
    private final Logger logger = new Logger("Workflow");

    public Workflow()
    {
        /*
        In order to skip steps without checking their hash, development mode has to be active.
         */

        steps.add(new tfprio.InitStaticVariables());
        /*
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

        steps.add(new lib.Tepic.Tepic());
        if (TFPRIO.configs.tepic.randomizeTfGeneMatrix.isSet() && TFPRIO.configs.tepic.randomizeTfGeneMatrix.get())
        {
            steps.add(new lib.Tepic.Randomize());
        }
        steps.add(new lib.Tepic.Postprocessing());
        steps.add(new lib.Plots.OpenRegionsViolinPlots());
         */

        if (TFPRIO.configs.tgene.pathToExecutable.isSet())
        {
            if (!TFPRIO.configs.mixOptions.mutuallyExclusive.get())
            {
                steps.add(new lib.Tgene.CreateGroups());
            }

            /*
            com2pose_lib.filter_target_genes_tgen();

            if (options_intern.tgen_self_regulatory)
            {
                com2pose_lib.integrate_self_regulatory_tgen();
            } */
        }
    }

    public boolean simulationSuccessful()
    {
        ExecutionTimeMeasurement timer = new ExecutionTimeMeasurement();
        logger.info("Start simulation.");
        boolean allWorked = true;
        for (ExecutableStep step : steps)
        {
            allWorked = allWorked && step.simulate();
        }
        logger.info("Finished simulation. Execution took " + timer.stopAndGetDeltaSeconds() + " seconds.");
        return allWorked;
    }

    public void run()
    {
        ExecutionTimeMeasurement timer = new ExecutionTimeMeasurement();
        logger.info("Start running executable steps.");
        for (ExecutableStep step : steps)
        {
            step.run();
        }
        logger.info(
                "Finished running executable steps. Execution took " + timer.stopAndGetDeltaSeconds() + " seconds.");
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
