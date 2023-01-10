package org.exbio.tfprio.steps.peakFiles;

import org.apache.commons.lang3.tuple.Pair;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.File;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.function.Function;
import java.util.stream.Collectors;

import static org.exbio.pipejar.util.FileManagement.extend;
import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

/**
 * HINT tool calculates footprints in ATAQ-seq data
 * Input:
 * - .narrowPeaks processed by MACS2
 * - .bam files
 * - genome (e.g. hg38)
 * - information about sequencing type (paired-, single-end)
 * Output:
 * - .bed footprint files
 */
public class HINT extends ExecutableStep {
    public final Map<String, Map<String, Collection<OutputFile>>> outputFiles = new HashMap<>();
    private final String exec = "rgt-hint";
    private final String mode = "footprinting";
    private final RequiredConfig<String> seqType = new RequiredConfig<>(Configs.mixOptions.seqType);
    private final RequiredConfig<File> bamDirectory = new RequiredConfig<>(Configs.hint.bam_directory);
    private final RequiredConfig<Boolean> paired = new RequiredConfig<>(Configs.hint.paired);
    private final RequiredConfig<String> genome = new RequiredConfig<>(Configs.hint.genome);
    private final Map<OutputFile, Pair<InputFile, InputFile>> bridge = new HashMap<>();
    private final Function<File, String> trimFile = (file) -> {
        String name = file.getName();
        return name.substring(0, name.lastIndexOf('.'));
    };

    /**
     * Init HINT with peak files
     *
     * @param peakFiles: Output of preprocessing steps for peakFiles
     */
    public HINT(Map<String, Map<String, Collection<OutputFile>>> peakFiles) {
        super(false, peakFiles.values().stream().flatMap(
                stringCollectionMap -> stringCollectionMap.values().stream()).flatMap(Collection::stream).collect(
                Collectors.toSet()));

        peakFiles.forEach((group, hmMap) -> {
            OutputFile d_groupIn = new OutputFile(inputDirectory, group);
            OutputFile d_groupOut = new OutputFile(outputDirectory, group);
            hmMap.forEach((hm, sampleFiles) -> {
                OutputFile d_hmIn = new OutputFile(d_groupIn, hm);
                OutputFile d_hmOut = new OutputFile(d_groupOut, hm);

                sampleFiles.forEach(bedFile -> {
                    String sampleName = trimFile.apply(bedFile);

                    OutputFile bamFile =
                            new OutputFile(extend(bamDirectory.get(), group, sampleName + ".bam").getAbsolutePath());
                    OutputFile bamIndex = new OutputFile(bamFile.getAbsolutePath() + ".bai");

                    InputFile inputBed = addInput(d_hmIn, bedFile);
                    InputFile inputBam = addInput(d_hmIn, bamFile);
                    if (bamIndex.exists()) {
                        addInput(d_hmIn, bamIndex);
                    }

                    OutputFile outputFile = addOutput(d_hmOut, inputBed.getName());

                    outputFiles.computeIfAbsent(group, k -> new HashMap<>()).computeIfAbsent(hm,
                            k -> new HashSet<>()).add(outputFile);
                    bridge.put(outputFile, Pair.of(inputBed, inputBam));
                });
            });
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((out, inputs) -> add(() -> {
                InputFile bed = inputs.getLeft();
                InputFile bam = inputs.getRight();

                // set sequencing type (one of: atac-seq; dnase-seq)
                String type = "--" + seqType.get();
                // basic hint call for given sequencing type
                ArrayList<String> command_args = new ArrayList<>(List.of(exec, mode, type));
                if (seqType.get().equals("dnase-seq")) {
                    // set genome
                    command_args.add("--bias-correction");
                }
                // paired or not
                if (paired.get()) {
                    command_args.add("--paired-end");
                }
                // organism
                command_args.add("--organism" + "=" + genome.get());
                // output location
                command_args.add("--output-location" + "=" + out.getParentFile().getAbsolutePath());
                // output prefix -> final name will be: ${prefix}.bed
                command_args.add("--output-prefix" + "=" + trimFile.apply(out));
                // add bam file
                command_args.add(bam.getAbsolutePath());
                // add peaks file
                command_args.add(bed.getAbsolutePath());
                executeAndWait(command_args, true);
                return true;
            }));
        }};
    }
}