package org.exbio.tfprio.steps.TEPIC;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.OptionalConfig;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.UsageConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.File;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import static org.exbio.pipejar.util.FileManagement.extend;
import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class TEPIC extends ExecutableStep {
    private final Map<InputFile, OutputFile> bridge = new HashMap<>();

    private final Map<String, Map<String, Collection<InputFile>>> chipSeqFiles = new HashMap<>();

    private final RequiredConfig<File> f_referenceGenome = new RequiredConfig<>(Configs.tepic.inputReferenceGenome);
    private final RequiredConfig<File> f_pwms = new RequiredConfig<>(Configs.tepic.PWMs);
    private final OptionalConfig<File> bedChromatinSignal =
            new OptionalConfig<>(Configs.tepic.bedChromatinSignal, false);
    private final OptionalConfig<Integer> columnBedfile = new OptionalConfig<>(Configs.tepic.columnBedfile, false);
    private final OptionalConfig<File> geneAnnotationFile =
            new OptionalConfig<>(Configs.tepic.geneAnnotationFile, false);
    private final OptionalConfig<Integer> windowSize = new OptionalConfig<>(Configs.tepic.windowSize, false);
    private final OptionalConfig<File> onlyDNasePeaks = new OptionalConfig<>(Configs.tepic.onlyDNasePeaks, false);
    private final OptionalConfig<Boolean> exponentialDecay =
            new OptionalConfig<>(Configs.tepic.exponentialDecay, false);
    private final OptionalConfig<Boolean> doNotNormalizePeakLength =
            new OptionalConfig<>(Configs.tepic.doNotNormalizePeakLength, false);
    private final OptionalConfig<Boolean> originalDecay = new OptionalConfig<>(Configs.tepic.originalDecay, false);
    private final OptionalConfig<File> psemsLengthFile = new OptionalConfig<>(Configs.tepic.psemsLengthFile, false);
    private final OptionalConfig<Boolean> entireGeneBody = new OptionalConfig<>(Configs.tepic.entireGeneBody, false);
    private final OptionalConfig<Boolean> doZip = new OptionalConfig<>(Configs.tepic.doZip, false);
    private final OptionalConfig<File> twoBitFile = new OptionalConfig<>(Configs.tepic.twoBitFile, false);
    private final OptionalConfig<Double> pValue = new OptionalConfig<>(Configs.tepic.pValue, false);
    private final OptionalConfig<Integer> maxMinutesPerChromosome =
            new OptionalConfig<>(Configs.tepic.maxMinutesPerChromosome, false);
    private final OptionalConfig<Boolean> chromosomePrefix =
            new OptionalConfig<>(Configs.tepic.chromosomePrefix, false);
    private final OptionalConfig<Boolean> transcriptBased = new OptionalConfig<>(Configs.tepic.transcriptBased, false);
    private final OptionalConfig<File> loopListFile = new OptionalConfig<>(Configs.tepic.loopListFile, false);
    private final OptionalConfig<Integer> loopWindows = new OptionalConfig<>(Configs.tepic.loopWindows, false);
    private final OptionalConfig<Boolean> onlyPeakFeatures =
            new OptionalConfig<>(Configs.tepic.onlyPeakFeatures, false);

    private final File executable;

    public TEPIC(Map<String, Map<String, Collection<OutputFile>>> latestChipSeq, OutputFile tepicDirectory) {
        super(false,
                latestChipSeq.values().stream().flatMap(x -> x.values().stream()).flatMap(Collection::stream).collect(
                        Collectors.toSet()), tepicDirectory);

        addInput(tepicDirectory);
        executable = extend(tepicDirectory, "Code", "TEPIC.sh");

        latestChipSeq.forEach((group, hmMap) -> {
            OutputFile dGroupOut = new OutputFile(outputDirectory, group);
            chipSeqFiles.put(group, new HashMap<>());
            OutputFile dGroupIn = new OutputFile(inputDirectory, group);
            hmMap.forEach((hm, samples) -> {
                OutputFile dHmOut = new OutputFile(dGroupOut, hm);
                chipSeqFiles.get(group).put(hm, new HashSet<>());
                OutputFile dHmIn = new OutputFile(dGroupIn, hm);
                samples.forEach(sample -> {
                    OutputFile sampleDir =
                            addOutput(dHmOut, sample.getName().substring(0, sample.getName().lastIndexOf('.')));
                    InputFile inputFile = addInput(dHmIn, sample);
                    bridge.put(inputFile, sampleDir);
                    chipSeqFiles.get(group).get(hm).add(inputFile);
                });
            });
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        String order = "gbopScdnawfmervkiqhs";
        return new HashSet<>() {{
            chipSeqFiles.forEach((group, hmMap) -> hmMap.forEach((hm, samples) -> {
                samples.forEach(sample -> add(() -> {
                    Map<Character, String> stringConfigs = new HashMap<>() {{
                        put('b', sample.getAbsolutePath());
                        put('o', bridge.get(sample).getAbsolutePath() + "/");
                        put('S', new File(bridge.get(sample), "trap_sequences.csv").getAbsolutePath());
                    }};

                    Map<Character, UsageConfig<?>> otherConfigs = new HashMap<>() {{
                        put('g', f_referenceGenome);
                        put('p', f_pwms);
                        put('d', bedChromatinSignal);
                        put('n', columnBedfile);
                        put('a', geneAnnotationFile);
                        put('w', windowSize);
                        put('f', onlyDNasePeaks);
                        put('e', exponentialDecay);
                        put('l', doNotNormalizePeakLength);
                        put('x', originalDecay);
                        put('m', psemsLengthFile);
                        put('y', entireGeneBody);
                        put('z', doZip);
                        put('r', twoBitFile);
                        put('v', pValue);
                        put('i', maxMinutesPerChromosome);
                        put('j', chromosomePrefix);
                        put('t', transcriptBased);
                        put('h', loopListFile);
                        put('s', loopWindows);
                        put('q', onlyPeakFeatures);
                    }};

                    otherConfigs.forEach((key, config) -> {
                        if (config.isSet()) {
                            if (config.get().getClass().equals(Boolean.class)) {
                                if (List.of('e', 'q').contains(key)) {
                                    // Boolean parameters
                                    stringConfigs.put(key, config.toString().toUpperCase());
                                } else {
                                    // Flags
                                    if ((boolean) config.get()) {
                                        stringConfigs.put(key, "");
                                    }
                                }
                            } else {
                                stringConfigs.put(key, config.get().toString());
                            }
                        }
                    });
                    String command = executable.getAbsolutePath() + " " +
                            order.chars().mapToObj(c -> (char) c).filter(stringConfigs::containsKey).map(c -> {
                                if (stringConfigs.containsKey(c)) {
                                    return "-" + c + " " + stringConfigs.get(c);
                                } else {
                                    return "";
                                }
                            }).collect(Collectors.joining(" "));

                    executeAndWait(command, true);

                    return true;
                }));
            }));
        }};
    }

    @Override
    protected boolean mayBeSkipped() {
        return false;
    }
}