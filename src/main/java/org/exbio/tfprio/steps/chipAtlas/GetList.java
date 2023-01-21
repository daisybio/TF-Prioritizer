package org.exbio.tfprio.steps.chipAtlas;

import org.apache.commons.io.IOUtils;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.File;
import java.net.URL;
import java.util.Collection;
import java.util.HashSet;
import java.util.concurrent.Callable;

public class GetList extends ExecutableStep<Configs> {
    public final OutputFile outputFile = addOutput("chip_atlas_file_list.csv");
    private final RequiredConfig<String> listUrl = new RequiredConfig<>(configs.chipAtlas.listUrl);

    public GetList(Configs configs) {
        super(configs);
    }

    @Override
    protected boolean doCreateFiles() {
        return false;
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                File zip = new File(outputDirectory, "chip_atlas_file_list.zip");

                int attempt = 1;

                while (!zip.exists()) {
                    logger.debug("Attempt " + attempt + " to download chip atlas list");
                    try {
                        IOUtils.copy(new URL(listUrl.get()), zip);
                    } catch (Exception e) {
                        attempt++;
                        if (attempt > 3) {
                            throw e;
                        }
                    }
                }
                logger.debug("Downloaded chip atlas list");

                String command =
                        "unzip " + zip.getAbsolutePath() + " -d " + outputFile.getParentFile().getAbsolutePath();

                Process unzip = Runtime.getRuntime().exec(command);
                int ret = unzip.waitFor();

                if (ret != 0) {
                    throw new RuntimeException("Unzip failed with return code " + ret);
                }

                return true;
            });
        }};
    }
}
