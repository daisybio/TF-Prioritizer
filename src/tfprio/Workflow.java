package tfprio;

import lib.ExecutableStep;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.ExecutionTimeMeasurement;
import util.Logger;

import java.io.File;
import java.util.*;

import static util.FileManagement.getFirstExisting;

public class Workflow
{
    private final List<ExecutableStep> steps = new ArrayList<>();
    private final Set<ExecutableStep> notExecutedSteps = new HashSet<>();
    private final Logger logger = new Logger("Workflow");

    public Workflow()
    {
        if (TFPRIO.configs.mixOptions.mutuallyExclusive.get() && !TFPRIO.configs.mixOptions.level.isSet())
        {
            TFPRIO.configs.mixOptions.level.setValue("SAMPLE_LEVEL");
        }

        steps.add(new tfprio.InitStaticVariables());
        steps.add(new lib.CheckChromosomes());

        if (TFPRIO.configs.mixOptions.level.isSet())
        {
            steps.add(new lib.MixOptions.SampleLevelPreprocess());
            steps.add(new lib.MixOptions.SampleLevel());

            if (TFPRIO.configs.mixOptions.level.get().equals("HM_LEVEL"))
            {
                steps.add(new lib.MixOptions.HmLevelPreprocess());
                steps.add(new lib.MixOptions.HmLevel());
            } else
            {
                for (ExecutableStep notExecutedStep : Arrays.asList(new lib.MixOptions.HmLevelPreprocess(),
                        new lib.MixOptions.HmLevel()))
                {
                    notExecutedSteps.add(notExecutedStep);
                    notExecutedStep.setNoExecutionReason(
                            TFPRIO.configs.mixOptions.level.getName() + " is not set to HM_LEVEL");
                }
            }
        } else
        {
            for (ExecutableStep notExecutedStep : Arrays.asList(new lib.MixOptions.SampleLevelPreprocess(),
                    new lib.MixOptions.SampleLevel(), new lib.MixOptions.HmLevelPreprocess(),
                    new lib.MixOptions.HmLevel()))
            {
                notExecutedSteps.add(notExecutedStep);
                notExecutedStep.setNoExecutionReason(TFPRIO.configs.mixOptions.level.getName() + " is not set");
            }
        }

        if (TFPRIO.configs.tepic.tfBindingSiteSearch.isSet() &&
                (TFPRIO.configs.tepic.tfBindingSiteSearch.get().equals("BETWEEN") ||
                        TFPRIO.configs.tepic.tfBindingSiteSearch.get().equals("EXCL_BETWEEN")))
        {
            steps.add(new lib.FootprintIntervals.CreateFootprintsBetweenPeaks());
        } else
        {
            ExecutableStep notExecutedStep = new lib.FootprintIntervals.CreateFootprintsBetweenPeaks();
            notExecutedStep.setNoExecutionReason(TFPRIO.configs.tepic.tfBindingSiteSearch.getName() + " does not " +
                    "equal BETWEEN or EXCL_BETWEEN");
            notExecutedSteps.add(notExecutedStep);
        }

        if (TFPRIO.configs.blacklist.bedFilePath.isSet())
        {
            steps.add(new lib.Blacklist.Blacklist());
        } else
        {
            ExecutableStep notExecutedStep = new lib.Blacklist.Blacklist();
            notExecutedStep.setNoExecutionReason(TFPRIO.configs.blacklist.bedFilePath.getName() + " is not set");
            notExecutedSteps.add(notExecutedStep);
        }

        if (TFPRIO.configs.mixOptions.mutuallyExclusive.get())
        {
            steps.add(new lib.MixOptions.MixMutuallyExclusive());
        } else
        {
            ExecutableStep notExecutedStep = new lib.MixOptions.MixMutuallyExclusive();
            notExecutedStep.setNoExecutionReason(
                    TFPRIO.configs.mixOptions.mutuallyExclusive.getName() + " is set to false");
            notExecutedSteps.add(notExecutedStep);
        }

        {
            steps.add(new lib.Deseq2.CreateDeseq2Scripts());
            steps.add(new lib.Deseq2.CreateTpmMappings());
            steps.add(new lib.Deseq2.CreateGenePositions());

            if (TFPRIO.configs.deSeq2.tpmFilter.isSet())
            {
                steps.add(new lib.Deseq2.PreprocessTpm());
            } else
            {
                ExecutableStep notExecutedStep = new lib.Deseq2.PreprocessTpm();
                notExecutedStep.setNoExecutionReason(TFPRIO.configs.deSeq2.tpmFilter.getName() + " is not set");
                notExecutedSteps.add(notExecutedStep);
            }
            steps.add(new lib.Deseq2.Deseq2());
            steps.add(new lib.Deseq2.PostProcessing());
        }

        steps.add(new lib.Tgene.Preprocess());
        steps.add(new lib.Tgene.RunTgene());
        steps.add(new lib.Tgene.Postprocessing());

        steps.add(new lib.Tepic.Tepic());
        if (TFPRIO.configs.tepic.randomizeTfGeneMatrix.isSet() && TFPRIO.configs.tepic.randomizeTfGeneMatrix.get())
        {
            steps.add(new lib.Tepic.Randomize());
        } else
        {
            ExecutableStep notExecutedStep = new lib.Tepic.Randomize();
            notExecutedStep.setNoExecutionReason(
                    TFPRIO.configs.tepic.randomizeTfGeneMatrix.getName() + " is set to false or " +
                            TFPRIO.configs.tepic.randomizeTfGeneMatrix.getName() + " is not set");
            notExecutedSteps.add(notExecutedStep);
        }
        steps.add(new lib.Tepic.Postprocessing());
        steps.add(new lib.Plots.OpenRegionsViolinPlots());

        if (TFPRIO.configs.tgene.pathToExecutable.isSet())
        {
            if (!TFPRIO.configs.mixOptions.mutuallyExclusive.get())
            {
                steps.add(new lib.Tgene.CreateGroups());
            } else
            {
                ExecutableStep notExecutedStep = new lib.Tgene.CreateGroups();
                notExecutedStep.setNoExecutionReason(
                        TFPRIO.configs.mixOptions.mutuallyExclusive.getName() + " is set to true");
                notExecutedSteps.add(notExecutedStep);
            }
            steps.add(new lib.Tgene.FilterTargetGenes());

            if (TFPRIO.configs.tgene.selfRegulatory.get())
            {
                steps.add(new lib.Tgene.SelfRegulatory());
            } else
            {
                ExecutableStep notExecutedStep = new lib.Tgene.SelfRegulatory();
                notExecutedStep.setNoExecutionReason(
                        TFPRIO.configs.tgene.selfRegulatory.getName() + " is set to false");
                notExecutedSteps.add(notExecutedStep);
            }
        } else
        {
            for (ExecutableStep notExecutedStep : Arrays.asList(new lib.Tgene.CreateGroups(),
                    new lib.Tgene.FilterTargetGenes(), new lib.Tgene.SelfRegulatory()))
            {
                notExecutedStep.setNoExecutionReason(TFPRIO.configs.tgene.pathToExecutable.getName() + " is not set");
                notExecutedSteps.add(notExecutedStep);
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

        steps.add(new lib.Logos.TfBindingLogoBiophysicalSequence());
        steps.add(new lib.Logos.PredictedBindingSites());

        if (TFPRIO.configs.chipAtlas.isEnabled.get())
        {
            steps.add(new lib.ChiPAtlas.GetDataList());
            steps.add(new lib.ChiPAtlas.GetData());
        } else
        {
            for (ExecutableStep notExecutedStep : Arrays.asList(new lib.ChiPAtlas.GetDataList(),
                    new lib.ChiPAtlas.GetData()))
            {
                notExecutedStep.setNoExecutionReason(TFPRIO.configs.chipAtlas.isEnabled.getName() + " is set to false");
                notExecutedSteps.add(notExecutedStep);
            }
        }

        steps.add(new lib.Igv.DcgTargetGenes());

        if (TFPRIO.configs.igv.importantLociAllPrioTf.isSet() && TFPRIO.configs.igv.importantLociAllPrioTf.isSet())
        {
            steps.add(new lib.Igv.ImportantLoci());
        } else
        {
            ExecutableStep notExecutedStep = new lib.Igv.ImportantLoci();
            notExecutedStep.setNoExecutionReason(
                    TFPRIO.configs.igv.importantLociAllPrioTf.getName() + " is not set or " +
                            TFPRIO.configs.igv.importantLociAllPrioTf.getName() + " is not set");
            notExecutedSteps.add(notExecutedStep);
        }
        if (TFPRIO.configs.igv.topLog2fc.isSet())
        {
            steps.add(new lib.Igv.TopLog2FC());
        } else
        {
            ExecutableStep notExecutedStep = new lib.Igv.TopLog2FC();
            notExecutedStep.setNoExecutionReason(TFPRIO.configs.igv.topLog2fc.getName() + " is not set");
            notExecutedSteps.add(notExecutedStep);
        }
        steps.add(new lib.DistributionAnalysis.CoOccurrenceAnalysis());
        steps.add(new lib.Report.Report());

        markNotGeneratedFileStructures();
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

    private void markNotGeneratedFileStructures()
    {
        for (ExecutableStep notExecutedStep : notExecutedSteps)
        {
            for (GeneratedFileStructure outputFileStructure : notExecutedStep.getCreatedFileStructure())
            {
                try
                {
                    outputFileStructure.deleteAndSetNoGenerationReason(notExecutedStep.getNoExecutionReason());
                } catch (ClassCastException e)
                {
                    logger.warn("File structure could not be marked as not created: " + outputFileStructure.getName());
                }
            }
        }
    }

    public static AbstractConfig<File> getLatestInputDirectory()
    {
        List<AbstractConfig<File>> inputDirectories =
                Arrays.asList(TFPRIO.configs.mixOptions.fileStructure.d_mutuallyExclusive_input,
                        TFPRIO.configs.blacklist.fileStructure.d_newInput,
                        TFPRIO.configs.mixOptions.fileStructure.d_footprintsBetweenPeaks,
                        TFPRIO.configs.mixOptions.fileStructure.d_hmMix,
                        TFPRIO.configs.mixOptions.fileStructure.d_sampleMix,
                        TFPRIO.configs.mixOptions.fileStructure.d_preprocessingCheckChr,
                        TFPRIO.configs.tepic.inputDirectory);

        return getFirstExisting(inputDirectories);
    }
}
