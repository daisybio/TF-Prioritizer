package org.exbio.tfprio;

import org.apache.commons.cli.ParseException;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.pipejar.pipeline.Workflow;
import org.exbio.tfprio.configs.Configs;
import org.exbio.tfprio.steps.Dynamite.*;
import org.exbio.tfprio.steps.EnsgSymbol;
import org.exbio.tfprio.steps.TEPIC.*;
import org.exbio.tfprio.steps.chipAtlas.CoOccurrenceAnalysis;
import org.exbio.tfprio.steps.chipAtlas.GetData;
import org.exbio.tfprio.steps.chipAtlas.GetList;
import org.exbio.tfprio.steps.distributionAnalysis.*;
import org.exbio.tfprio.steps.igv.DistributionTargetGenes;
import org.exbio.tfprio.steps.logos.*;
import org.exbio.tfprio.steps.peakFiles.*;
import org.exbio.tfprio.steps.plots.*;
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
        ExecutableStep.setAcceptAllInputs();

        new TF_Prioritizer(args);
    }

    @Override
    protected Configs createConfigs() {
        return new Configs();
    }

    protected void buildFlow() {
        // read peaks input
        InitPeaks initPeaks = add(new InitPeaks());
        Map<String, Map<String, Collection<OutputFile>>> latestPeakFiles = initPeaks.outputFiles;

        // perform HINT preprocessing for ATAC-/DNASE-seq input data
        if (List.of("atac-seq", "dnase-seq").contains(Configs.mixOptions.seqType.get())) {
            HINT hint = add(new HINT(latestPeakFiles));
            latestPeakFiles = hint.outputFiles;
        }

        CheckChromosomes checkChromosomes = add(new CheckChromosomes(latestPeakFiles));
        latestPeakFiles = checkChromosomes.outputFiles;

        if (Configs.mixOptions.perform.get()) {
            MixSamples mixSamples = add(new MixSamples(latestPeakFiles));
            latestPeakFiles = mixSamples.outputFiles;
        }

        if (List.of("BETWEEN", "EXCL_BETWEEN").contains(Configs.mixOptions.tfBindingSiteSearch.get())) {
            CreateFootprintsBetweenPeaks createFootprintsBetweenPeaks =
                    add(new CreateFootprintsBetweenPeaks(latestPeakFiles));
            latestPeakFiles = createFootprintsBetweenPeaks.outputFiles;
        }

        if (Configs.mixOptions.blackListPath.isSet()) {
            Blacklist blacklist = add(new Blacklist(latestPeakFiles));
            latestPeakFiles = blacklist.outputFiles;
        }


        GetChromosomeLengths getChromosomeLengths = add(new GetChromosomeLengths());

        if (Configs.mixOptions.mixMutuallyExclusive.get()) {
            MixMutuallyExclusive mixMutuallyExclusive = add(new MixMutuallyExclusive(latestPeakFiles));
            latestPeakFiles = mixMutuallyExclusive.outputFiles;
        }

        FilterEnsgs filterEnsgs = add(new FilterEnsgs());

        MergeCounts mergeCounts = add(new MergeCounts(filterEnsgs.outputFile));

        FetchGeneInfo fetchGeneInfo = add(new FetchGeneInfo(filterEnsgs.cleanFile));

        EnsgSymbol ensgSymbol = add(new EnsgSymbol(filterEnsgs.cleanFile));

        CreateBatchFile createBatchFile = add(new CreateBatchFile(mergeCounts.outputFiles));
        MeanExpression meanCounts = add(new MeanExpression(mergeCounts.outputFiles));
        CreatePairings createPairings = add(new CreatePairings(mergeCounts.outputFiles));

        CalculateTPM calculateTPM = add(new CalculateTPM(meanCounts.outputFiles, fetchGeneInfo.outputFile));

        Uplift uplift = add(new Uplift(fetchGeneInfo.outputFile));
        FilterENdb filterENdb = add(new FilterENdb());
        UpliftENdb upliftENdb = add(new UpliftENdb(filterENdb.outputFile));

        if (Configs.deSeq2.tpmFilter.isSet()) {
            FilterExpression filterExpression = add(new FilterExpression(calculateTPM.outputFiles));
        }

        DeSeq2 deSeq2 = add(new DeSeq2(createPairings.outputFiles, createBatchFile.outputFile));
        DeSeqPostprocessing deSeqPostprocessing = add(new DeSeqPostprocessing(deSeq2.outputFiles));

        Map<String, Map<String, OutputFile>> tgeneFiles = new HashMap<>();

        if (Configs.tGene.executable.isSet()) {
            TGenePreprocess tGenePreprocess = add(new TGenePreprocess());
            TGeneExtractRegions tGeneExtractRegions = add(new TGeneExtractRegions(tGenePreprocess.outputFile));
            TGene tGene = add(new TGene(latestPeakFiles, tGenePreprocess.outputFile));
            TGenePostprocessing tGenePostprocessing =
                    add(new TGenePostprocessing(meanCounts.outputFiles, tGene.outputFiles));
            tgeneFiles = tGenePostprocessing.outputFiles;
        }
        PreprocessReferenceGenome preprocessReferenceGenome = add(new PreprocessReferenceGenome());
        TEPIC tepic = add(new TEPIC(latestPeakFiles, preprocessReferenceGenome.outputFile));
        Map<String, Map<String, Collection<OutputFile>>> tepicFiles = tepic.outputFiles;

        if (Configs.tepic.randomize.isSet() && Configs.tepic.randomize.get()) {
            TEPICRandomize tepicRandomize = add(new TEPICRandomize(tepic.outputFiles));
            tepicFiles = tepicRandomize.outputFiles;
        }

        ExtractAffinities extractAffinities = add(new ExtractAffinities(tepicFiles));
        ExtractSequences extractSequences = add(new ExtractSequences(tepicFiles));

        CreateBindingRegionsBedFiles createBindingRegionsBedFiles =
                add(new CreateBindingRegionsBedFiles(extractSequences.outputFiles));
        CalculateMeanAffinities calculateMeanAffinities =
                add(new CalculateMeanAffinities(extractAffinities.outputFiles));
        CalculateAffinityRatios calculateAffinityRatios =
                add(new CalculateAffinityRatios(calculateMeanAffinities.outputFiles));
        OpenRegionsViolinPlots openRegionsViolinPlots = add(new OpenRegionsViolinPlots(tepicFiles));
        FindAnalyzableTFs findAnalyzableTFs = add(new FindAnalyzableTFs(calculateMeanAffinities.outputFiles));

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
        ExtractRegressionCoefficients extractRegressionCoefficients =
                add(new ExtractRegressionCoefficients(runDynamite.outputFiles));

        FilterRegressionCoefficients filterRegressionCoefficients =
                add(new FilterRegressionCoefficients(extractRegressionCoefficients.outputFiles));
        ThresholdPlots thresholdPlots = add(new ThresholdPlots(filterRegressionCoefficients.outputFiles));

        GroupStages groupStages = add(new GroupStages(filterRegressionCoefficients.outputFiles));
        PlotGroupedStages plotGroupedStages = add(new PlotGroupedStages(groupStages.outputFiles));
        AnalyzeGroupLevel analyzeGroupLevel =
                add(new AnalyzeGroupLevel(calculateTPM.outputFiles, meanCounts.outputFiles, groupStages.outputFiles,
                        ensgSymbol.outputFile, findAnalyzableTFs.outputFile));
        AnalyzeHmLevel analyzeHmLevel = add(new AnalyzeHmLevel(analyzeGroupLevel.outputFiles));

        org.exbio.tfprio.steps.plots.TopKTargetGenes topKTargetGenes =
                add(new org.exbio.tfprio.steps.plots.TopKTargetGenes(groupStages.outputFiles,
                        calculateMeanAffinities.outputFiles));

        Preprocessing daPreprocessing = add(new Preprocessing(groupStages.outputFiles));
        RunDistributionAnalysis runDistributionAnalysis =
                add(new RunDistributionAnalysis(daPreprocessing.outputFile, ensgSymbol.outputFile,
                        meanCounts.outputFiles, deSeqPostprocessing.outputFiles,
                        extractRegressionCoefficients.outputFiles, calculateMeanAffinities.outputFiles));
        CreateBackground createBackground = add(new CreateBackground(runDistributionAnalysis.outputFiles));
        CreatePlots createPlots =
                add(new CreatePlots(runDistributionAnalysis.outputFiles, createBackground.outputFiles));
        ExtractStatRank extractStatRank = add(new ExtractStatRank(createPlots.statFiles));
        CalculateDcgRank calculateDcgRank = add(new CalculateDcgRank(extractStatRank.outputFiles));

        org.exbio.tfprio.steps.distributionAnalysis.TopKTargetGenes daTopKTargetGenes =
                add(new org.exbio.tfprio.steps.distributionAnalysis.TopKTargetGenes(calculateDcgRank.outputFile,
                        calculateMeanAffinities.outputFiles));
        CreateHeatmaps createHeatmaps = add(new CreateHeatmaps(deSeq2.normalizedCounts, daTopKTargetGenes.outputFiles,
                createBatchFile.outputFile, ensgSymbol.outputFile));
        BiophysicalModel biophysicalModel = add(new BiophysicalModel(calculateDcgRank.outputFile));
        BiophysicalLogo biophysicalLogo = add(new BiophysicalLogo(biophysicalModel.outputFile));

        Jaspar jaspar = add(new Jaspar(calculateDcgRank.outputFile));
        ExtractBindingSites extractBindingSites =
                add(new ExtractBindingSites(calculateDcgRank.outputFile, extractSequences.outputFiles));
        FilterBindingSites filterBindingSites = add(new FilterBindingSites(extractBindingSites.outputFiles));
        PredictedLogo predictedLogo = add(new PredictedLogo(filterBindingSites.outputFiles));

        OutputFile chipAtlasDirectory = null;
        if (Configs.chipAtlas.enabled.isSet() && Configs.chipAtlas.enabled.get()) {
            GetList getList = add(new GetList());
            GetData getData = add(new GetData(getList.outputFile, calculateDcgRank.outputFile));
            chipAtlasDirectory = getData.outputFile;

            CoOccurrenceAnalysis coOccurrenceAnalysis = add(new CoOccurrenceAnalysis(chipAtlasDirectory));
        }

        DistributionTargetGenes distributionTargetGenes =
                add(new DistributionTargetGenes(ensgSymbol.outputFile, fetchGeneInfo.outputFile,
                        calculateDcgRank.outputFile, chipAtlasDirectory, daTopKTargetGenes.outputFiles, tepicFiles,
                        createBindingRegionsBedFiles.outputFiles));
        distributionTargetGenes.setUnderDevelopment();
    }
}
