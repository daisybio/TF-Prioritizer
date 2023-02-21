package org.exbio.tfprio;

import org.apache.commons.cli.ParseException;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.Workflow;
import org.exbio.tfprio.configs.Configs;
import org.exbio.tfprio.steps.Dynamite.*;
import org.exbio.tfprio.steps.EnsgSymbol;
import org.exbio.tfprio.steps.TEPIC.*;
import org.exbio.tfprio.steps.chipAtlas.*;
import org.exbio.tfprio.steps.distributionAnalysis.*;
import org.exbio.tfprio.steps.igv.DistributionTargetGenes;
import org.exbio.tfprio.steps.igv.ImportantLoci;
import org.exbio.tfprio.steps.igv.TopLog2fc;
import org.exbio.tfprio.steps.logos.*;
import org.exbio.tfprio.steps.peakFiles.*;
import org.exbio.tfprio.steps.plots.GroupStages;
import org.exbio.tfprio.steps.plots.OpenRegionsViolinPlots;
import org.exbio.tfprio.steps.report.ReportCompression;
import org.exbio.tfprio.steps.report.ReportCreation;
import org.exbio.tfprio.steps.report.ReportPreprocessing;
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
        // read peaks input
        InitPeaks initPeaks = add(new InitPeaks(configs));
        Map<String, Map<String, Collection<OutputFile>>> latestPeakFiles = initPeaks.outputFiles;

        // perform HINT preprocessing for ATAC-/DNASE-seq input data
        if (List.of("atac-seq", "dnase-seq").contains(configs.inputConfigs.seqType.get())) {
            HINT hint = add(new HINT(configs, latestPeakFiles));
            latestPeakFiles = hint.outputFiles;
        }


        CheckChromosomes checkChromosomes = add(new CheckChromosomes(configs, latestPeakFiles));
        latestPeakFiles = checkChromosomes.outputFiles;

        if (configs.mixOptions.perform.get()) {
            MixSamples mixSamples = add(new MixSamples(configs, latestPeakFiles));
            latestPeakFiles = mixSamples.outputFiles;
        }

        if (List.of("BETWEEN", "EXCL_BETWEEN").contains(configs.mixOptions.tfBindingSiteSearch.get())) {
            CreateFootprintsBetweenPeaks createFootprintsBetweenPeaks =
                    add(new CreateFootprintsBetweenPeaks(configs, latestPeakFiles));
            latestPeakFiles = createFootprintsBetweenPeaks.outputFiles;
        }

        if (configs.mixOptions.blackListPath.isSet()) {
            Blacklist blacklist = add(new Blacklist(configs, latestPeakFiles));
            latestPeakFiles = blacklist.outputFiles;
        }


        if (configs.mixOptions.mixMutuallyExclusive.get()) {
            MixMutuallyExclusive mixMutuallyExclusive = add(new MixMutuallyExclusive(configs, latestPeakFiles));
            latestPeakFiles = mixMutuallyExclusive.outputFiles;
        }

        FilterEnsgs filterEnsgs = add(new FilterEnsgs(configs));

        Map<String, OutputFile> rawCounts;

        if (configs.inputConfigs.rnaSeq.get().isDirectory()) {
            MergeCounts mergeCounts = add(new MergeCounts(configs, filterEnsgs.outputFile));
            rawCounts = mergeCounts.outputFiles;
        } else {
            SplitCounts splitCounts = add(new SplitCounts(configs, filterEnsgs.outputFile));
            rawCounts = splitCounts.outputFiles;
        }

        FetchGeneInfo fetchGeneInfo = add(new FetchGeneInfo(configs, filterEnsgs.cleanFile));

        EnsgSymbol ensgSymbol = add(new EnsgSymbol(configs, filterEnsgs.cleanFile));

        CreateBatchFile createBatchFile = add(new CreateBatchFile(configs, rawCounts));
        MeanExpression meanCounts = add(new MeanExpression(configs, rawCounts));
        CreatePairings createPairings = add(new CreatePairings(configs, rawCounts));

        CalculateTPM calculateTPM = add(new CalculateTPM(configs, meanCounts.outputFiles, fetchGeneInfo.outputFile));

        Uplift uplift = add(new Uplift(configs, fetchGeneInfo.outputFile));
        FilterENdb filterENdb = add(new FilterENdb(configs));
        UpliftENdb upliftENdb = add(new UpliftENdb(configs, filterENdb.outputFile));

        if (configs.deSeq2.tpmFilter.isSet()) {
            FilterExpression filterExpression = add(new FilterExpression(configs, calculateTPM.outputFiles));
        }

        DeSeq2 deSeq2 = add(new DeSeq2(configs, createPairings.outputFiles, createBatchFile.outputFile));
        DeSeqPostprocessing deSeqPostprocessing = add(new DeSeqPostprocessing(configs, deSeq2.outputFiles));

        Map<String, Map<String, OutputFile>> tgeneFiles = new HashMap<>();

        if (configs.tGene.executable.isSet()) {
            TGenePreprocess tGenePreprocess = add(new TGenePreprocess(configs));
            TGeneExtractRegions tGeneExtractRegions = add(new TGeneExtractRegions(configs, tGenePreprocess.outputFile));
            TGene tGene = add(new TGene(configs, latestPeakFiles, tGenePreprocess.outputFile));
            TGenePostprocessing tGenePostprocessing =
                    add(new TGenePostprocessing(configs, meanCounts.outputFiles, tGene.outputFiles));
            tgeneFiles = tGenePostprocessing.outputFiles;
        }
        PreprocessReferenceGenome preprocessReferenceGenome = add(new PreprocessReferenceGenome(configs));
        TEPIC tepic = add(new TEPIC(configs, latestPeakFiles, preprocessReferenceGenome.outputFile));
        Map<String, Map<String, Collection<OutputFile>>> tepicFiles = tepic.outputFiles;

        if (configs.tepic.randomize.isSet() && configs.tepic.randomize.get()) {
            TEPICRandomize tepicRandomize = add(new TEPICRandomize(configs, tepic.outputFiles));
            tepicFiles = tepicRandomize.outputFiles;
        }

        ExtractAffinities extractAffinities = add(new ExtractAffinities(configs, tepicFiles));
        ExtractSequences extractSequences = add(new ExtractSequences(configs, tepicFiles));

        CreateBindingRegionsBedFiles createBindingRegionsBedFiles =
                add(new CreateBindingRegionsBedFiles(configs, extractSequences.outputFiles));
        CalculateMeanAffinities calculateMeanAffinities =
                add(new CalculateMeanAffinities(configs, extractAffinities.outputFiles));
        CalculateAffinityRatios calculateAffinityRatios =
                add(new CalculateAffinityRatios(configs, calculateMeanAffinities.outputFiles));
        OpenRegionsViolinPlots openRegionsViolinPlots = add(new OpenRegionsViolinPlots(configs, tepicFiles));

        if (configs.tGene.executable.isSet()) {
            if (!configs.mixOptions.mixMutuallyExclusive.get()) {
                add(new TGeneCreateGroups());
            }

            if (configs.tGene.selfRegulatory.isSet() && configs.tGene.selfRegulatory.get()) {
                add(new TGeneSelfRegulatory());
            }
        }

        IntegrateData integrateData =
                add(new IntegrateData(configs, calculateAffinityRatios.outputFiles, deSeqPostprocessing.getOutputs()));
        PrepareForClassification prepareForClassification =
                add(new PrepareForClassification(configs, integrateData.outputFiles));
        RunDynamite runDynamite = add(new RunDynamite(configs, prepareForClassification.outputFiles));
        ExtractRegressionCoefficients extractRegressionCoefficients =
                add(new ExtractRegressionCoefficients(configs, runDynamite.outputFiles));
        FilterRegressionCoefficients filterRegressionCoefficients =
                add(new FilterRegressionCoefficients(configs, extractRegressionCoefficients.outputFiles));

        GroupStages groupStages = add(new GroupStages(configs, filterRegressionCoefficients.outputFiles));

        Preprocessing daPreprocessing = add(new Preprocessing(configs, groupStages.outputFiles));
        RunDistributionAnalysis runDistributionAnalysis =
                add(new RunDistributionAnalysis(configs, daPreprocessing.outputFile, ensgSymbol.outputFile,
                        meanCounts.outputFiles, deSeqPostprocessing.outputFiles,
                        filterRegressionCoefficients.outputFiles, calculateMeanAffinities.outputFiles));
        CreateBackground createBackground = add(new CreateBackground(configs, runDistributionAnalysis.outputFiles));
        CreatePlots createPlots =
                add(new CreatePlots(configs, runDistributionAnalysis.outputFiles, createBackground.outputFiles));
        ExtractStatRank extractStatRank = add(new ExtractStatRank(configs, createPlots.statFiles));
        CalculateDcgPerHm calculateDcgPerHm = add(new CalculateDcgPerHm(configs, extractStatRank.outputFiles));
        CalculateDcgRank calculateDcgRank = add(new CalculateDcgRank(configs, extractStatRank.outputFiles));

        org.exbio.tfprio.steps.distributionAnalysis.TopKTargetGenes daTopKTargetGenes =
                add(new org.exbio.tfprio.steps.distributionAnalysis.TopKTargetGenes(configs,
                        calculateDcgRank.outputFile, calculateMeanAffinities.outputFiles));
        CreateHeatmaps createHeatmaps =
                add(new CreateHeatmaps(configs, deSeq2.normalizedCounts, daTopKTargetGenes.outputFiles,
                        createBatchFile.outputFile, ensgSymbol.outputFile));
        BiophysicalModel biophysicalModel = add(new BiophysicalModel(configs, calculateDcgRank.outputFile));
        BiophysicalLogo biophysicalLogo = add(new BiophysicalLogo(configs, biophysicalModel.outputFile));

        Jaspar jaspar = add(new Jaspar(configs, calculateDcgRank.outputFile));
        ExtractBindingSites extractBindingSites =
                add(new ExtractBindingSites(configs, calculateDcgRank.outputFile, extractSequences.outputFiles));
        FilterBindingSites filterBindingSites = add(new FilterBindingSites(configs, extractBindingSites.outputFiles));
        PredictedLogo predictedLogo = add(new PredictedLogo(configs, filterBindingSites.outputFiles));

        OutputFile chipAtlasDirectory = null;
        OutputFile coOccurrence = null;
        OutputFile confusionMatrixesDir = null;
        if (configs.chipAtlas.enabled.isSet() && configs.chipAtlas.enabled.get()) {
            GetList getList = add(new GetList(configs));
            GetData getData = add(new GetData(configs, getList.outputFile, calculateDcgRank.outputFile));
            chipAtlasDirectory = getData.outputFile;

            CoOccurrenceAnalysis coOccurrenceAnalysis = add(new CoOccurrenceAnalysis(configs, chipAtlasDirectory));
            coOccurrence = coOccurrenceAnalysis.outputFile;

            CoOccurrenceBindingEnergies coOccurrenceBindingEnergies =
                    add(new CoOccurrenceBindingEnergies(configs, extractSequences.outputFiles,
                            coOccurrenceAnalysis.merged));
            GetChromosomeLengths getChromosomeLengths = add(new GetChromosomeLengths(configs));
            ConfusionMatrixes confusionMatrixes =
                    add(new ConfusionMatrixes(configs, getData.outputFile, getChromosomeLengths.outputFile,
                            createBindingRegionsBedFiles.outputFiles));
            confusionMatrixesDir = confusionMatrixes.outputFile;
        }

        DistributionTargetGenes distributionTargetGenes =
                add(new DistributionTargetGenes(configs, ensgSymbol.outputFile, fetchGeneInfo.outputFile,
                        calculateDcgRank.outputFile, chipAtlasDirectory, daTopKTargetGenes.outputFiles, tepicFiles,
                        createBindingRegionsBedFiles.outputFiles));

        Map<String, OutputFile> importantLociFiles;
        if (configs.igv.importantLoci.isSet()) {
            ImportantLoci importantLoci =
                    add(new ImportantLoci(configs, fetchGeneInfo.outputFile, ensgSymbol.outputFile,
                            chipAtlasDirectory));
            importantLociFiles = importantLoci.outputFiles;
        } else {
            importantLociFiles = new HashMap<>();
        }

        TopLog2fc topLog2fc =
                add(new TopLog2fc(configs, deSeq2.outputFiles, fetchGeneInfo.outputFile, ensgSymbol.outputFile,
                        chipAtlasDirectory));

        ReportPreprocessing reportPreprocessing =
                add(new ReportPreprocessing(configs, ensgSymbol.outputFile, calculateDcgPerHm.outputFiles,
                        deSeq2.outputFiles, meanCounts.outputFiles, calculateTPM.outputFiles,
                        biophysicalLogo.outputFile, jaspar.outputFile, createHeatmaps.outputFiles,
                        distributionTargetGenes.outputFiles, createPlots.plotFiles,
                        extractRegressionCoefficients.outputFiles, coOccurrence, importantLociFiles,
                        topLog2fc.outputFiles, confusionMatrixesDir));

        ReportCreation reportCreation = add(new ReportCreation(configs, reportPreprocessing.outputFile));
        ReportCompression reportCompression = add(new ReportCompression(configs, reportCreation.outputFile));
    }
}
