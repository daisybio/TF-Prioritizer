package tfprio;

import lib.ExecutableStep;

import java.util.ArrayList;
import java.util.List;

public class Workflow
{
    private final List<ExecutableStep> steps = new ArrayList<>();

    public Workflow()
    {
        steps.add(new lib.CheckChromosomes());
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

        steps.add(new lib.Deseq2.CreateDeseq2Scripts());
        steps.add(new lib.Deseq2.CreateTpmMappings());
        steps.add(new lib.Deseq2.CreateGenePositions());

        if (TFPRIO.configs.deSeq2.tpmFilter.isSet())
        {
            steps.add(new lib.Deseq2.PreprocessTpm());
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

    public void run()
    {
        for (ExecutableStep step : steps)
        {
            step.run();
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
