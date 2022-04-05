package tfprio;

import lib.ExecutableStep;
import util.ExecutionTimeMeasurement;
import util.Logger;

import java.util.*;

public class Workflow
{
    private final List<ExecutableStep> steps = new ArrayList<>();
    private final Logger logger = new Logger("Workflow");

    public Workflow()
    {
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
            steps.add(new lib.Tgene.Preprocess());
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

        if (TFPRIO.configs.tgene.pathToExecutable.isSet())
        {
            if (!TFPRIO.configs.mixOptions.mutuallyExclusive.get())
            {
                steps.add(new lib.Tgene.CreateGroups());
            }
            steps.add(new lib.Tgene.filterTargetGenes());

            if (TFPRIO.configs.tgene.selfRegulatory.get())
            {
                steps.add(new lib.Tgene.SelfRegulatory());
            }
        }

        steps.add(new lib.Dynamite.Preprocessing());
        steps.add(new lib.Dynamite.RunDynamite());

        steps.add(new lib.Plots.GroupPlots());
        steps.add(new lib.Plots.AnalyzeData());
        steps.add(new lib.Plots.TopKTargetGenes());

        steps.add(new lib.DistributionAnalysis.Preprocessing());
        steps.add(new lib.DistributionAnalysis.RunDistributionAnalysis());
        steps.add(new lib.DistributionAnalysis.CreatePlots());
        steps.add(new lib.DistributionAnalysis.CalculateDcgRank());
        steps.add(new lib.DistributionAnalysis.getTopKTargetGenes());
        steps.add(new lib.DistributionAnalysis.GenerateTargetGenesHeatmaps());

        steps.add(new lib.Logos.TfBindingLogoBiophysicalSequence());*/
        steps.add(new lib.Logos.PredictedBindingSites());
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
            System.gc(); // Trigger garbage collection
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
