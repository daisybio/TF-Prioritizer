package org.exbio.tfprio;

import org.exbio.pipejar.configs.ConfigModuleCollection;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.pipejar.pipeline.ExecutionManager;
import org.exbio.pipejar.steps.ConcatenateFiles;
import org.exbio.tfprio.configs.Configs;
import org.exbio.tfprio.steps.chipSeq.*;
import org.exbio.tfprio.steps.rnaSeq.*;
import org.exbio.tfprio.steps.tGene.TGene;
import org.exbio.tfprio.steps.tGene.TGeneExtractRegions;
import org.exbio.tfprio.steps.tGene.TGenePostprocessing;
import org.exbio.tfprio.steps.tGene.TGenePreprocess;
import org.exbio.tfprio.util.ArgParser;

import java.io.File;
import java.io.IOException;
import java.util.*;

import static org.exbio.pipejar.util.FileManagement.extend;

public class TF_Prioritizer {
    private static final Collection<ExecutableStep> steps = new HashSet<>();
    /**
     * The {@link ConfigModuleCollection} object which is referenced by all the {@link ExecutableStep} instances
     */
    public static Configs configs;
    public static File workingDirectory;

    public static void main(String[] args) throws Exception {
        ArgParser argParser = new ArgParser(args);
        init(argParser);
        buildFlow();
        execute();
    }

    private static void init(ArgParser argParser) throws IOException {
        File configFile = argParser.getConfigFile();
        workingDirectory = argParser.getWorkingDirectory();
        ExecutionManager.workingDirectory = new OutputFile(extend(workingDirectory, "output").getAbsolutePath());
        ExecutionManager.setThreadNumber(5);

        configs = new Configs();

        if (!configs.merge(configFile) || !configs.validate()) {
            System.exit(1);
        }
        configs.save(extend(workingDirectory, "configs.json"));
    }

    protected static void buildFlow() {
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

        GetChromosomeLengths getChromosomeLengths = add(new GetChromosomeLengths());

        if (Configs.mixOptions.mixMutuallyExclusive.get()) {
            MixMutuallyExclusive mixMutuallyExclusive = add(new MixMutuallyExclusive(latestChipSeq));
            latestChipSeq = mixMutuallyExclusive.outputFiles;
        }

        FetchGeneInfo fetchGeneInfo = add(new FetchGeneInfo());
        ConcatenateFiles concatenateGeneInfo = add(new ConcatenateFiles(fetchGeneInfo.getOutputs()));
        MergeCounts mergeCounts = add(new MergeCounts());

        CreateBatchFile createBatchFile = add(new CreateBatchFile(mergeCounts.outputFiles));
        CalculateTPM calculateTPM = add(new CalculateTPM(mergeCounts.outputFiles, concatenateGeneInfo.outputFile));
        MeanExpression meanCounts = add(new MeanExpression(mergeCounts.outputFiles));
        CreatePairings createPairings = add(new CreatePairings(mergeCounts.outputFiles));

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

        //TEPIC tepic = add(new TEPIC(latestChipSeq));
    }

    private static void execute() {
        ExecutionManager manager = new ExecutionManager(steps);
        manager.run();
    }

    private static <T extends ExecutableStep> T add(T step) {
        steps.add(step);
        return step;
    }
}
