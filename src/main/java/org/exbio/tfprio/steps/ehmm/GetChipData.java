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
import java.util.concurrent.ForkJoinPool;
import java.util.stream.Collectors;

import static org.exbio.pipejar.util.FileManagement.makeSureFileExists;
import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class GetChipData extends ExecutableStep<Configs> {
    public final OutputFile outputFile = addOutput("background.bed");
    public final OutputFile bedDir = addOutput("bed");
    public final OutputFile bamDir = addOutput("bam");
    private final InputFile chipFileList;
    private final InputFile chromosomeLengths;
    private final List<String> tissues = new RequiredConfig<>(configs.chipAtlas.tissueTypes).get();
    private final String genomeVersion = new RequiredConfig<>(configs.inputConfigs.genome).get();
    private final Set<String> antigenClassKeys= new RequiredConfig<>(configs.ehmm.antigenClassKeys).get();
    private final String threshold = new RequiredConfig<>(configs.ehmm.threshold).get();

    public GetChipData(Configs configs, OutputFile chipFileList, OutputFile chromosomeLengths){
        super(configs, false, chipFileList);
        this.chipFileList = addInput(chipFileList);
        this.chromosomeLengths = addInput(chromosomeLengths);
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
        public final String key;

        public ChipListEntry(String line){
            String[] array = line.split(",");
            this.fileName = array[0];
            this.genomeAssembly = array[1];
            this.antigenClass = array[2];
            this.antigen = array[3];
            this.cellTypeClass = array[4];
            this.cellType = array[5];
            this.threshold = array[6];
            this.fileUrl = array[array.length-1];
            this.key = this.antigen.equals("") ? this.antigenClass: this.antigenClass + "_" + this.antigen;
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
                logger.debug("filtering for params: genome: {}, antigenClassKeys: {}, tissues: {}, threshold: {}",
                        genomeVersion, antigenClassKeys, tissues, threshold);
                Set<ChipListEntry> validEntries = Files.lines(Path.of(chipFileList.getAbsolutePath()))
                        .skip(1)
                        .map(ChipListEntry::new)
                        .filter(entry -> entry.genomeAssembly.equalsIgnoreCase(genomeVersion))
                        .filter(entry -> antigenClassKeys.stream()
                                .anyMatch(antigenClass -> antigenClass.equalsIgnoreCase(entry.key)))
                        .filter(entry -> tissues.stream().
                                anyMatch(tissue -> tissue.equalsIgnoreCase(entry.cellTypeClass)))
                        .filter(entry -> Objects.equals(entry.threshold, threshold))
                        .filter(entry -> entry.cellType.equals(""))
                        .collect(Collectors.toSet());
                int nEntries = validEntries.size();
                if (nEntries > 0) {
                    logger.debug("Found {} matching bed files in ChipAtlas: {}",
                            nEntries, validEntries.stream().map(ChipListEntry::name).collect(Collectors.toSet()));
                } else {
                    throw new RuntimeException("No matching bed files found in ChipAtlas for given filters");
                }
                ForkJoinPool threadPool = new ForkJoinPool(nEntries);
                try {
                    threadPool.submit(() ->
                            validEntries.parallelStream().forEach(entry -> {
                                File bedFile = addOutput(bedDir, entry.name().replace(" ", "_"));
                                int attempt = 1;
                                while (!bedFile.exists()) {
                                    logger.debug("Downloading chip atlas bed file: " + bedFile.getName() + " (Attempt: " + attempt + ")");
                                    try {
                                        makeSureFileExists(bedFile);
                                        IOUtils.copy(new URL(entry.fileUrl), bedFile);

                                        long size = Files.size(bedFile.toPath());

                                        if (size < 300) {
                                            logger.warn("File {} is too small ({} bytes)", bedFile, size);
                                            throw new IOException("File too small");
                                        }
                                    } catch (IOException e) {
                                        attempt++;
                                        if (attempt > 3) {
                                            throw new RuntimeException("Could not download file " + entry.fileUrl);
                                        }
                                    }
                                }
                                logger.trace("Simplifying ChipAtlas file {}", bedFile.getAbsolutePath());
                                OutputFile simpleBed = new OutputFile(bedDir, "tmp.bed");
                                // simplify data from chipAtlas, causes errors in bam creation
                                try (BufferedReader br = new BufferedReader(new FileReader(bedFile));
                                     BufferedWriter bw = new BufferedWriter(new FileWriter(simpleBed))) {
                                    br.readLine();  // skip header
                                    for (String line = br.readLine(); line != null; line = br.readLine()) {
                                        String[] tsvData = line.split("\t");
                                        String chr = tsvData[0];
                                        String auxData = tsvData[3];
                                        // remove chr prefix
                                        if (chr.startsWith("chr")) chr = chr.substring("chr".length());
                                        // simplify aux data
                                        auxData = auxData.substring(0, auxData.indexOf(";"));
                                        tsvData[0] = chr;
                                        tsvData[3] = auxData;
                                        bw.write(String.join("\t", tsvData));
                                        bw.newLine();
                                    }
                                } catch (IOException e) {
                                    throw new RuntimeException("ChipAtlas bed file does not exist or cannot be opened.");
                                }
                                // replace complex bed with simplified version
                                try {
                                    Files.move(simpleBed.toPath(), bedFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
                                } catch (IOException e) {
                                    throw new RuntimeException("Failed to replace " + bedFile.getAbsolutePath() +
                                            " with " + simpleBed.getAbsolutePath());
                                }
                                logger.trace("Building associated bam file for {}", bedFile.getAbsolutePath());
                                OutputFile bamFile = addOutput(bamDir, FilenameUtils.getBaseName(bedFile.getPath()) + ".bam");
                                String cmd = String.join(" ",
                                        "bedToBam",
                                        "-i", bedFile.getAbsolutePath(),
                                        "-g", chromosomeLengths.getAbsolutePath(),
                                        "|", "samtools", "sort", "-",
                                        ">", bamFile.getAbsolutePath(),
                                        ";", "samtools", "index", bamFile.getAbsolutePath());
                                try {
                                    executeAndWait(List.of("/bin/sh", "-c", cmd), false);
                                } catch (IOException e) {
                                    throw new RuntimeException("BedToBam creation failed");
                                }
                                try {
                                    Files.write(Paths.get(outputFile.getAbsolutePath()),
                                            (Iterable<String>) Files.lines(Path.of(bedFile.getAbsolutePath()))::iterator,
                                            StandardOpenOption.APPEND);
                                } catch (IOException e) {
                                    throw new RuntimeException("Appending to merged bed file failed");
                                }
                            }));
                } finally {
                    threadPool.shutdown();
                }
                return true;
            });
        }};
    }
}
