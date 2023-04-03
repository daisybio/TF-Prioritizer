package org.exbio.tfprio.steps.ehmm;

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.*;
import java.net.URL;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.exbio.pipejar.util.FileManagement.deleteFileStructure;
import static org.exbio.pipejar.util.FileManagement.makeSureFileExists;
import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class GetChipData extends ExecutableStep<Configs> {
    public final OutputFile outputFile = addOutput("background.bed");
    public final OutputFile bedDir = addOutput("bed");
    public final OutputFile bamDir = addOutput("bam");
    private final InputFile chipFileList;
    private final InputFile chromosomeLengths;
    private final RequiredConfig<List<String>> tissues = new RequiredConfig<>(configs.chipAtlas.tissueTypes);
    private final RequiredConfig<String> genomeVersion = new RequiredConfig<>(configs.inputConfigs.genome);
    private String seqType = new RequiredConfig<>(configs.inputConfigs.seqType).get();
    private final Set<String> histoneModifications;

    public GetChipData(Configs configs, OutputFile chipFileList,
                       Map<String, Map<String, Collection<OutputFile>>> peakFiles, OutputFile chromosomeLengths){
        super(configs, false, chipFileList);
        this.chipFileList = addInput(chipFileList);
        this.chromosomeLengths = addInput(chromosomeLengths);

        // no specific modifications or antigens to filter ChipAtlas for
        if (!this.seqType.equals("chip-seq")) {
            this.histoneModifications = Set.of("ALL");
        } else {
            this.seqType = "histone";
            this.histoneModifications = peakFiles.values().stream().findFirst().orElseThrow().keySet();
        }
    }

    private static class ChipListEntry {
        public final String fileName;
        public final String genomeAssembly;
        public final String antigenClass;
        public final String antigen;
        public final String cellTypeClass;
        public final String cellType;
        public final String threshold;
        public final String fileUrl;

        private String chooseAntigenClass(String antigen){
            return antigen.equalsIgnoreCase("Histone") ? "chip-seq": antigen.toLowerCase();
        }

        public ChipListEntry(String line){
            String[] array = line.split(",");
            this.fileName = array[0];
            this.genomeAssembly = array[1];
            this.antigenClass = this.chooseAntigenClass(array[2]);
            this.antigen = array[3];
            this.cellTypeClass = array[4];
            this.cellType = array[5];
            this.threshold = array[6];
            this.fileUrl = array[array.length-1];
        }

        public String name(){
            return String.join("_",
                    List.of((this.antigen.equals("") ? "ALL": this.antigen),
                            this.threshold,
                            this.cellTypeClass,
                            this.antigenClass)) + ".bed";
        }
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                logger.trace("Starting to read the ChipAtlas file list");
                Stream<ChipListEntry> validEntries = Files.lines(Path.of(chipFileList.getAbsolutePath()))
                        .skip(1)
                        .map(ChipListEntry::new)
                        .filter(entry -> entry.genomeAssembly.equalsIgnoreCase(genomeVersion.get()))
                        .filter(entry -> entry.antigenClass.equalsIgnoreCase(seqType))
                        .filter(entry -> tissues.get().stream().
                                anyMatch(tissue -> tissue.equalsIgnoreCase(entry.cellTypeClass)))
                        .filter(entry -> histoneModifications.stream()
                                .anyMatch(hm -> hm.equalsIgnoreCase(entry.antigen)));
                int nEntries = validEntries.collect(Collectors.toSet()).size();
                if (nEntries > 0) {
                    logger.debug("Found {} matching bed files in ChipAtlas", nEntries);
                } else {
                    throw new RuntimeException("No matching bed files found in ChipAtlas for given filters");
                }
                validEntries.forEach(entry -> add(() -> {
                    File bedFile = addOutput(bedDir, entry.name());
                    int attempt = 1;
                    while (!bedFile.exists()) {
                        logger.debug("Attempt " + attempt + " to download chip atlas bed file: " + bedFile.getName());
                        try {
                            makeSureFileExists(bedFile);
                            IOUtils.copy(new URL(entry.fileUrl), bedFile);

                            long size = Files.size(bedFile.toPath());

                            if (size < 300) {
                                logger.warn("File {} is too small ({} bytes)", bedFile, size);
                                throw new IOException("File too small");
                            }
                        } catch (IOException e) {
                            deleteFileStructure(bedFile);
                            attempt++;
                            if (attempt > 3) {
                                throw new RuntimeException("Could not download file " + entry.fileUrl);
                            }
                        }
                    }
                    OutputFile simpleBed = new OutputFile(bedDir, "tmp.bed");
                    // simplify aux data column from chipAtlas, causes errors in bam creation
                    try (BufferedReader br = new BufferedReader(new FileReader(bedFile));
                         BufferedWriter bw = new BufferedWriter(new FileWriter(simpleBed))) {
                        for (String line = br.readLine(); line != null; line = br.readLine()) {
                            bw.write(line.substring(line.indexOf(";")));
                            bw.newLine();
                        }
                    } catch (IOException e) {
                        throw new RuntimeException("ChipAtlas bed file does not exist or cannot be opened.");
                    }
                    // replace complex bed with simplified version
                    Files.move(simpleBed.toPath(), bedFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
                    String base = FilenameUtils.getBaseName(bedFile.getPath());
                    // create associated bam file
                    OutputFile bamFile = addOutput(bamDir, base +".bam");
                    String bedToBam = String.join(" ",
                            "bedToBam",
                            "-i", bedFile.getAbsolutePath(),
                            "-g", chromosomeLengths.getAbsolutePath(),
                            ">", bamFile.getAbsolutePath());
                    try {
                        executeAndWait(bedToBam, false);
                    } catch (IOException e) {
                        throw new RuntimeException("Failed to convert bed to bam:\n" + e.getMessage());
                    }
                    // sorted bam file
                    OutputFile sortedBamFile = addOutput(bamDir, base + ".sorted.bam");
                    String sortBam = String.join(" ",
                            "samtools", "sort", bamFile.getAbsolutePath(),
                            ">", sortedBamFile.getAbsolutePath());
                    // replace bam with sorted file
                    Files.move(sortedBamFile.toPath(), bamFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
                    try {
                        executeAndWait(sortBam, false);
                    } catch (IOException e) {
                        throw new RuntimeException("Failed to sort bam file:\n" + e.getMessage());
                    }
                    // index bam file
                    String indexBam = String.join(" ",
                            "samtools", "index", bamFile.getAbsolutePath());
                    try {
                        executeAndWait(indexBam, false);
                    } catch (IOException e) {
                        throw new RuntimeException("Failed to index bam file:\n" + e.getMessage());
                    }
                    // merge to one bed file
                    Files.write(Paths.get(outputFile.getAbsolutePath()),
                            (Iterable<String>) Files.lines(Path.of(bedFile.getAbsolutePath()))::iterator,
                            StandardOpenOption.APPEND);
                    return true;
                }));
                return true;
            });
        }};
    }
}
