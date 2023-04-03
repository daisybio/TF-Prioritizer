package org.exbio.tfprio.steps.ehmm;

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.IOException;
import java.net.URI;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.Collection;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class EPDNew extends ExecutableStep<Configs> {
    public final OutputFile bedFile = addOutput("promoters.bed");
    public final OutputFile bamFile = addOutput("promoters.bam");
    private final InputFile chromosomeLengths;
    private final URI downloadURI;

    public EPDNew(Configs configs, OutputFile chromosomeLengths){
        super(configs, false, chromosomeLengths);
        this.chromosomeLengths = addInput(chromosomeLengths);
        String genome = new RequiredConfig<>(configs.inputConfigs.genome).get();
        URI baseURI = URI.create(new RequiredConfig<>(configs.epdNew.baseURL).get());
        Map<String, String> genomeMap = new RequiredConfig<>(configs.epdNew.genomeMap).get();
        if (!genomeMap.containsKey(genome)) throw new RuntimeException("Genome version not recognized in EPDnew");
        this.downloadURI = baseURI.resolve(genomeMap.get(genome));
    }


    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                String promoterBed = IOUtils.toString(downloadURI.toURL().openStream(), StandardCharsets.UTF_8);
                Files.writeString(bedFile.toPath(), promoterBed.replace(" ", "\t"));
                String base = FilenameUtils.getBaseName(bedFile.getPath());
                // create associated bam file
                OutputFile bamFile = addOutput(base +".bam");
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
                OutputFile sortedBamFile = addOutput(base + ".sorted.bam");
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
                return true;
            });
        }};
    }
}
