package org.exbio.tfprio;

import org.exbio.pipejar.configs.ConfigModuleCollection;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.pipejar.pipeline.ExecutionManager;
import org.exbio.pipejar.steps.ConcatenateFiles;
import org.exbio.tfprio.configs.Configs;
import org.exbio.tfprio.steps.CheckAnnotationChromosomes;
import org.exbio.tfprio.steps.chipSeq.*;
import org.exbio.tfprio.steps.rnaSeq.*;
import org.exbio.tfprio.util.ArgParser;

import java.io.File;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import static org.exbio.pipejar.util.FileManagement.extend;

public class TFPRIO {
    /**
     * The {@link ConfigModuleCollection} object which is referenced by all the {@link ExecutableStep} instances
     */
    public static Configs configs;
    public static File workingDirectory;
    public static File sourceDirectory;

    public static void main(String[] args) throws Exception {
        ArgParser argParser = new ArgParser(args);
        execute(argParser);
    }

    protected static void execute(ArgParser argParser) throws Exception {
        File configFile = argParser.getConfigFile();
        workingDirectory = argParser.getWorkingDirectory();
        sourceDirectory = argParser.getSourceDirectory();
        ExecutionManager.workingDirectory = new OutputFile(extend(workingDirectory, "output").getAbsolutePath());
        ExecutionManager.setThreadNumber(5);

        configs = new Configs();

        if (!configs.merge(configFile) || !configs.validate()) {
            System.exit(1);
        }
        configs.save(extend(workingDirectory, "configs.json"));

        Collection<ExecutableStep> steps = new HashSet<>();

        CheckChromosomes checkChromosomes = new CheckChromosomes();
        steps.add(checkChromosomes);

        Map<String, Map<String, Collection<OutputFile>>> latestChipSeq = checkChromosomes.outputFiles;

        if (Configs.mixOptions.perform.get()) {
            MixSamples mixSamples = new MixSamples(checkChromosomes.outputFiles);
            steps.add(mixSamples);
            latestChipSeq = mixSamples.outputFiles;
        }

        if (List.of("BETWEEN", "EXCL_BETWEEN").contains(Configs.mixOptions.tfBindingSiteSearch.get())) {
            CreateFootprintsBetweenPeaks createFootprintsBetweenPeaks = new CreateFootprintsBetweenPeaks(latestChipSeq);
            steps.add(createFootprintsBetweenPeaks);
            latestChipSeq = createFootprintsBetweenPeaks.outputFiles;
        }

        if (Configs.mixOptions.blackListPath.isSet()) {
            Blacklist blacklist = new Blacklist(latestChipSeq);
            steps.add(blacklist);
            latestChipSeq = blacklist.outputFiles;
        }

        GetChromosomeLengths getChromosomeLengths = new GetChromosomeLengths();
        steps.add(getChromosomeLengths);

        if (Configs.mixOptions.mixMutuallyExclusive.get()) {
            MixMutuallyExclusive mixMutuallyExclusive = new MixMutuallyExclusive(latestChipSeq);
            steps.add(mixMutuallyExclusive);
            latestChipSeq = mixMutuallyExclusive.outputFiles;
        }


        MergeCounts mergeCounts = new MergeCounts();
        steps.add(mergeCounts);
        CreatePairings createPairings = new CreatePairings(mergeCounts.outputFiles);
        steps.add(createPairings);

        CreateBatchFile createBatchFile = new CreateBatchFile(mergeCounts.outputFiles);
        steps.add(createBatchFile);


        FetchGeneInfo fetchGeneInfo = new FetchGeneInfo();
        steps.add(fetchGeneInfo);

        ConcatenateFiles concatenateGeneInfo = new ConcatenateFiles(fetchGeneInfo.getOutputs());
        steps.add(concatenateGeneInfo);

        CalculateTPM calculateTPM = new CalculateTPM(createPairings.getOutputs(), concatenateGeneInfo.outputFile);
        steps.add(calculateTPM);
        Collection<OutputFile> latestRnaSeq = calculateTPM.getOutputs();

        Uplift uplift = new Uplift(concatenateGeneInfo.outputFile);
        steps.add(uplift);

        FilterENdb filterENdb = new FilterENdb();
        steps.add(filterENdb);

        UpliftENdb upliftENdb = new UpliftENdb(filterENdb.outputFile);
        steps.add(upliftENdb);

        if (Configs.deSeq2.tpmFilter.isSet()) {
            FilterTPM filterTPM = new FilterTPM(calculateTPM.getOutputs());
            steps.add(filterTPM);
            latestRnaSeq = filterTPM.getOutputs();
        }

        DeSeq2 deSeq2 = new DeSeq2(latestRnaSeq, createBatchFile.outputFile);
        steps.add(deSeq2);

        DeSeqPostprocessing deSeqPostprocessing = new DeSeqPostprocessing(deSeq2.outputFiles);
        steps.add(deSeqPostprocessing);

        CheckAnnotationChromosomes checkAnnotationChromosomes = new CheckAnnotationChromosomes();
        steps.add(checkAnnotationChromosomes);

        ExecutionManager manager = new ExecutionManager(steps);
        manager.run();
    }
}
