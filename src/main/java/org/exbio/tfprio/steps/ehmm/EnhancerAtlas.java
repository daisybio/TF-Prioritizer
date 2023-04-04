package org.exbio.tfprio.steps.ehmm;

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.*;
import java.net.URI;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;


public class EnhancerAtlas extends ExecutableStep<Configs> {
    public final OutputFile bedDir = addOutput("bed");
    public final OutputFile bamDir = addOutput("bam");
    public final OutputFile outputFile = addOutput("allEnhancers.bed");
    private final InputFile chromosomeLengths;
    private final String enhancerHTML = new RequiredConfig<>(configs.enhancerAtlas.enhancerHTMLString).get();
    private final URI enhancerBaseURI = URI.create(new RequiredConfig<>(configs.enhancerAtlas.enhancerBaseURLString).get());
    private final RequiredConfig<String> genome = new RequiredConfig<>(configs.inputConfigs.genome);
    private final RequiredConfig<Map<String, String>> groups = new RequiredConfig<>(configs.inputConfigs.sameStages);
    private final RequiredConfig<List<String>> tissues = new RequiredConfig<>(configs.chipAtlas.tissueTypes);

    private final RequiredConfig<Map<String, String>> genomeToKeyMap = new RequiredConfig<>(configs.enhancerAtlas.genomeToKeyMap);
    private final RequiredConfig<Map<String, String>> enhancerVersionMap = new RequiredConfig<>(configs.enhancerAtlas.enhancerVersionMap);
    private final InputFile script;

    public EnhancerAtlas(Configs configs, OutputFile chromosomeLengths) {
        super(configs, false, chromosomeLengths);
        this.chromosomeLengths = addInput(chromosomeLengths);
        script = addInput(getClass().getResourceAsStream("uplift.py"), "uplift.py");
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            // get HTML of enhancerAtlas download page
            String enhancerHTMLString;
            try {
                enhancerHTMLString = IOUtils.toString(new URL(enhancerHTML).openStream(), StandardCharsets.UTF_8);
            } catch (IOException e) {
                throw new RuntimeException("Failed to fetch EnhancerAtlas download page");
            }
            // search for links of organism, tissues, and groups
            String genomeKey = genomeToKeyMap.get().get(genome.get());
            Collection<String> searchKeys = tissues.get();
            searchKeys.addAll(groups.get().keySet());
            String tissueKeys = "(" + String.join("|", searchKeys) + ")";
            Pattern p = Pattern.compile("href.*" + "(" + genomeKey + ".?" + tissueKeys + ".*\\.bed).*onclick", Pattern.CASE_INSENSITIVE);
            Matcher m = p.matcher(enhancerHTMLString);
            List<String> bedLinks = new ArrayList<>();
            // add links that match the given parameters
            while(m.find()) bedLinks.add(m.group(1));
            if (bedLinks.size()==0){
                throw new RuntimeException("No matching enhancers found in EnhancerAtlas for pattern: " + p
                + "(List of available cell lines can be found here: " + enhancerHTML + ")");
            } else {
                logger.debug("Found {} matching files in EnhancerAtlas: {}", bedLinks.size(), bedLinks);
            }
            // download beds, uplift to given genome version, create bam files, and merge bed to one
            String enhancerVersion = enhancerVersionMap.get().get(genomeKey);
            boolean isSameGenome = enhancerVersion.equalsIgnoreCase(genome.get());
            bedLinks.forEach(bedLink -> add(() -> {
                String tissueFile = bedLink.split("/")[1];
                File bedFile = new File(bedDir, tissueFile);
                logger.debug("downloading file: {}", enhancerBaseURI.resolve(bedLink));
                IOUtils.copy(enhancerBaseURI.resolve(bedLink).toURL(), bedFile);
                if (!isSameGenome) {
                    logger.debug("uplifting enhancer bed file {}", bedFile);
                    File liftedBed = new File(bedDir, FilenameUtils.getBaseName(tissueFile) + "_lifted.bed");
                    try {
                        executeAndWait("python3 " + script + " " +
                                String.join(" ",
                                        bedFile.getAbsolutePath(),
                                        liftedBed.getAbsolutePath(),
                                        enhancerVersion,
                                        genome.get(), "0", "1", "2"), true);
                    } catch (IOException e) {
                        throw new RuntimeException("Failed to uplift bed file: " + bedFile + "\ninternal error:\n" + e);
                    }
                    bedFile = liftedBed;
                }
                File checkChromosomes = new File(bedDir, FilenameUtils.getBaseName(tissueFile) + "_noChr.bed");
                // simplify aux data column from chipAtlas and remove chr prefix, causes errors in bam creation
                try (BufferedReader br = new BufferedReader(new FileReader(bedFile));
                     BufferedWriter bw = new BufferedWriter(new FileWriter(checkChromosomes))) {
                    for (String line = br.readLine(); line != null; line = br.readLine()) {
                        String base = line.startsWith("chr") ? line.substring("chr".length()) : line;
                        bw.write(base);
                        bw.newLine();
                    }
                } catch (IOException e) {
                    throw new RuntimeException("Could not remove chr prefix from EnhancerAtlas bed file");
                }
                // replace bedFile
                Files.move(checkChromosomes.toPath(), bedFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
                // create bam files
                OutputFile bamFile = addOutput(bamDir, FilenameUtils.getBaseName(bedFile.getPath()) +".bam");
                String cmd = String.join(" ",
                        "bedToBam",
                        "-i", bedFile.getAbsolutePath(),
                        "-g", chromosomeLengths.getAbsolutePath(),
                        "|", "samtools", "sort", "-",
                        ">", bamFile.getAbsolutePath(),
                        ";", "samtools", "index", bamFile.getAbsolutePath());
                executeAndWait(List.of("/bin/sh", "-c", cmd), false);
                Files.write(Paths.get(outputFile.getAbsolutePath()),
                        (Iterable<String>) Files.lines(Path.of(bedFile.getAbsolutePath()))::iterator,
                        StandardOpenOption.APPEND);
                return true;
            }));
        }};
    }
}
