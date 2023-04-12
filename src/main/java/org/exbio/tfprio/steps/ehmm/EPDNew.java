package org.exbio.tfprio.steps.ehmm;

import org.apache.commons.io.IOUtils;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.net.URI;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Collection;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.Callable;

public class EPDNew extends ExecutableStep<Configs> {
    public final OutputFile bedFile = addOutput("promoters.bed");
    private final URI downloadURI;

    public EPDNew(Configs configs){
        super(configs);
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
                Files.writeString(bedFile.toPath(), promoterBed
                        .replace(" ", "\t")
                        .replace("chr", ""));
                return true;
            });
        }};
    }
}
