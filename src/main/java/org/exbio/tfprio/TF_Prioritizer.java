package org.exbio.tfprio;

import org.apache.commons.cli.ParseException;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.Workflow;
import org.exbio.pipejar.steps.ConcatenateFiles;
import org.exbio.tfprio.configs.Configs;
import org.exbio.tfprio.steps.Dynamite.FilterRegressionCoefficients;
import org.exbio.tfprio.steps.Dynamite.IntegrateData;
import org.exbio.tfprio.steps.Dynamite.PrepareForClassification;
import org.exbio.tfprio.steps.Dynamite.RunDynamite;
import org.exbio.tfprio.steps.EnsgSymbol;
import org.exbio.tfprio.steps.TEPIC.*;
import org.exbio.tfprio.steps.chipSeq.*;
import org.exbio.tfprio.steps.plots.GroupStages;
import org.exbio.tfprio.steps.plots.OpenRegionsViolinPlots;
import org.exbio.tfprio.steps.plots.PlotGroupedStages;
import org.exbio.tfprio.steps.plots.ThresholdPlots;
import org.exbio.tfprio.steps.rnaSeq.*;
import org.exbio.tfprio.steps.tGene.*;

import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class TF_Prioritizer extends Workflow<Configs> {
    public TF_Prioritizer(String[] args) throws IOException, ParseException {
        super(args);
    }

    public static void main(String[] args) throws Exception {
        new TF_Prioritizer(args);
    }

    @Override
    protected Configs createConfigs() {
        return new Configs();
    }

    protected void buildFlow() {
        CheckChromosomes checkChromosomes = add(new CheckChromosomes());
        Map<String, Map<String, Collection<OutputFile>>> latestChipSeq = checkChromosomes.outputFiles;

        if (Configs.mixOptions.perform.get()) {
            MixSamples mixSamples = add(new MixSamples(checkChromosomes.outputFiles));
            latestChipSeq = mixSamples.outputFiles;
        }

        if (List.of("BETWEEN", "EXCL_BETWEEN").contains(Configs.mixOptions.tfBindingSiteSearch.get())) {
            CreateFootprintsBetweenPeaks createFootprintsBetweenPeaks =
                    add(new CreateFootprintsBetweenPeaks(latestChipSeq));
            latestChipSeq = createFootprintsBetweenPeaks.outputFiles;
        }

        if (Configs.mixOptions.blackListPath.isSet()) {
            Blacklist blacklist = add(new Blacklist(latestChipSeq));
            latestChipSeq = blacklist.outputFiles;
        }

        EnsgSymbol ensgSymbol = add(new EnsgSymbol());

        GetChromosomeLengths getChromosomeLengths = add(new GetChromosomeLengths());

        if (Configs.mixOptions.mixMutuallyExclusive.get()) {
            MixMutuallyExclusive mixMutuallyExclusive = add(new MixMutuallyExclusive(latestChipSeq));
            latestChipSeq = mixMutuallyExclusive.outputFiles;
        }

        FetchGeneInfo fetchGeneInfo = add(new FetchGeneInfo());
        ConcatenateFiles concatenateGeneInfo = add(new ConcatenateFiles(fetchGeneInfo.getOutputs()));
        MergeCounts mergeCounts = add(new MergeCounts());

        CreateBatchFile createBatchFile = add(new CreateBatchFile(mergeCounts.outputFiles));
        MeanExpression meanCounts = add(new MeanExpression(mergeCounts.outputFiles));
        CreatePairings createPairings = add(new CreatePairings(mergeCounts.outputFiles));

        CalculateTPM calculateTPM = add(new CalculateTPM(meanCounts.outputFiles, concatenateGeneInfo.outputFile));

        Uplift uplift = add(new Uplift(concatenateGeneInfo.outputFile));
        FilterENdb filterENdb = add(new FilterENdb());
        UpliftENdb upliftENdb = add(new UpliftENdb(filterENdb.outputFile));

        if (Configs.deSeq2.tpmFilter.isSet()) {
            FilterExpression filterExpression = add(new FilterExpression(calculateTPM.outputFiles));
        }

        DeSeq2 deSeq2 = add(new DeSeq2(createPairings.getOutputs(), createBatchFile.outputFile));
        DeSeqPostprocessing deSeqPostprocessing = add(new DeSeqPostprocessing(deSeq2.outputFiles));

        Map<String, Map<String, OutputFile>> tgeneFiles = new HashMap<>();

        if (Configs.tGene.executable.isSet()) {
            TGenePreprocess tGenePreprocess = add(new TGenePreprocess());
            TGeneExtractRegions tGeneExtractRegions = add(new TGeneExtractRegions(tGenePreprocess.outputFile));
            TGene tGene = add(new TGene(latestChipSeq, tGenePreprocess.outputFile));
            TGenePostprocessing tGenePostprocessing =
                    add(new TGenePostprocessing(meanCounts.outputFiles, tGene.outputFiles));
            tgeneFiles = tGenePostprocessing.outputFiles;
        }
        TEPIC tepic = add(new TEPIC(latestChipSeq));
        Map<String, Map<String, Collection<OutputFile>>> tepicFiles = tepic.outputFiles;

        if (Configs.tepic.randomize.isSet() && Configs.tepic.randomize.get()) {
            TEPICRandomize tepicRandomize = add(new TEPICRandomize(tepic.outputFiles));
            tepicFiles = tepicRandomize.outputFiles;
        }

        CreateBindingRegionsBedFiles createBindingRegionsBedFiles = add(new CreateBindingRegionsBedFiles(tepicFiles));
        CalculateMeanAffinities calculateMeanAffinities = add(new CalculateMeanAffinities(tepicFiles));
        CalculateAffinityRatios calculateAffinityRatios =
                add(new CalculateAffinityRatios(calculateMeanAffinities.outputFiles));
        OpenRegionsViolinPlots openRegionsViolinPlots = add(new OpenRegionsViolinPlots(tepicFiles));

        if (Configs.tGene.executable.isSet()) {
            if (!Configs.mixOptions.mixMutuallyExclusive.get()) {
                add(new TGeneCreateGroups());
            }

            if (Configs.tGene.selfRegulatory.isSet() && Configs.tGene.selfRegulatory.get()) {
                add(new TGeneSelfRegulatory());
            }
        }

        IntegrateData integrateData =
                add(new IntegrateData(calculateAffinityRatios.outputFiles, deSeqPostprocessing.getOutputs()));
        PrepareForClassification prepareForClassification =
                add(new PrepareForClassification(integrateData.outputFiles));
        RunDynamite runDynamite = add(new RunDynamite(prepareForClassification.outputFiles));

        FilterRegressionCoefficients filterRegressionCoefficients =
                add(new FilterRegressionCoefficients(runDynamite.outputFiles));
        ThresholdPlots thresholdPlots = add(new ThresholdPlots(filterRegressionCoefficients.outputFiles));
        GroupStages groupStages = add(new GroupStages(filterRegressionCoefficients.outputFiles));
        PlotGroupedStages plotGroupedStages = add(new PlotGroupedStages(groupStages.outputFiles));
    }
}
