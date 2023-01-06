package org.exbio.tfprio.steps.peakFiles;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.File;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

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
    private final RequiredConfig<File> inputDirectory = new RequiredConfig<>(Configs.hint.bam_directory);
    private final RequiredConfig<Boolean> paired = new RequiredConfig<>(Configs.hint.paired);
    private final RequiredConfig<String> genome = new RequiredConfig<>(Configs.hint.genome);
    private final Map<String, File> bamFiles = new HashMap<>();
    private final Map<InputFile, OutputFile> bridge = new HashMap<>();

    /**
     * Init HINT with peak files
     *
     * @param peakFiles: Output of preprocessing steps for peakFiles
     */
    public HINT(Map<String, Map<String, Collection<OutputFile>>> peakFiles) {
        super(false, peakFiles.values().stream().flatMap(
                stringCollectionMap -> stringCollectionMap.values().stream()).flatMap(Collection::stream).collect(
                Collectors.toSet()));

        // manage bam files as input
        setBamFiles();
        // manage output files
        peakFiles.forEach((group, hmMap) -> {
            outputFiles.put(group, new HashMap<>());
            OutputFile d_groupOut = addOutput(group);
            hmMap.forEach((hm, sampleFiles) -> {
                outputFiles.get(group).put(hm, new HashSet<>());
                OutputFile d_hmOut = addOutput(d_groupOut, hm);

                sampleFiles.forEach(sampleFile -> {
                    InputFile inputFile = addInput(sampleFile);

                    OutputFile outputFile = addOutput(d_hmOut, inputFile.getName());
                    outputFiles.get(group).get(hm).add(outputFile);
                    bridge.put(inputFile, outputFile);
                });
            });
        });
    }

    @Override
    protected boolean mayBeSkipped() {
        return false;
    }

    private String trimFile(String file) {
        return file.substring(0, file.indexOf('.'));
    }

    /**
     * Add bam files as input
     */
    private void setBamFiles() {
        InputFile input = addInput(inputDirectory);
        // create different groups
        Arrays.stream(Objects.requireNonNull(input.listFiles())).forEach(d_group -> {
            // separation folder (atac-seq)
            Arrays.stream(Objects.requireNonNull(d_group.listFiles())).forEach(d_hm -> {
                // add bam files
                Arrays.stream(Objects.requireNonNull(d_hm.listFiles())).forEach(bam -> {
                    File bam_file = new File(bam.getAbsolutePath());
                    String sample = trimFile(bam_file.getName());
                    bamFiles.put(sample, bam_file);
                });
            });
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((inputFile, out) -> add(() -> {
                File bam = bamFiles.get(trimFile(inputFile.getName()));
                // set sequencing type (one of: atac-seq; dnase-seq)
                String type = "--" + seqType.get();
                // basic hint call for given sequencing type
                ArrayList<String> command_args = new ArrayList<>(List.of(exec, mode, type));
                // paired or not
                if (paired.get()) {
                    command_args.add("--paired-end");
                }
                // organism
                command_args.add("--organism" + "=" + genome.get());
                // output location
                command_args.add("--output-location" + "=" + out.getParentFile().getAbsolutePath());
                // output prefix -> final name will be: ${prefix}.bed
                command_args.add("--output-prefix" + "=" + trimFile(out.getName()));
                // add bam file
                command_args.add(bam.getAbsolutePath());
                // add peaks file
                command_args.add(inputFile.getAbsolutePath());
                executeAndWait(command_args, true);
                return true;
            }));
        }};
    }
}