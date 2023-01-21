package org.exbio.tfprio.steps.logos;

import org.apache.commons.io.IOUtils;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.Collection;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import static org.exbio.pipejar.util.FileManagement.readLines;

public class Jaspar extends ExecutableStep<Configs> {
    public final OutputFile outputFile = addOutput("jaspar");
    private final RequiredConfig<String> jasparURL = new RequiredConfig<>(configs.internalConfigs.jasparURL);
    private final RequiredConfig<String> jasparLogoUrl = new RequiredConfig<>(configs.internalConfigs.jasparLogoUrl);
    private final InputFile dcgFile;

    public Jaspar(Configs configs, OutputFile dcgFile) {
        super(configs, false, dcgFile);
        this.dcgFile = addInput(dcgFile);
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            Set<String> tfNames;
            try {
                tfNames = readLines(dcgFile).stream().skip(1).map(line -> line.split("\t")[0]).collect(
                        Collectors.toSet());
            } catch (IOException e) {
                throw new RuntimeException(e);
            }

            Map<String, Set<String>> tfMatrixIds;

            try (BufferedReader reader = new BufferedReader(
                    new InputStreamReader(new URL(jasparURL.get()).openStream()))) {
                tfMatrixIds = reader.lines().filter(line -> line.startsWith(">")).map(line -> line.substring(1)).map(
                        line -> line.split("\t")).collect(Collectors.groupingBy(line -> line[1].toUpperCase(),
                        Collectors.mapping(line -> line[0], Collectors.toSet())));
            } catch (IOException e) {
                throw new RuntimeException(e);
            }

            tfNames.retainAll(tfMatrixIds.keySet());

            tfNames.forEach(tf -> {
                File tfDirectory = new File(outputFile, tf);
                tfDirectory.mkdirs();

                tfMatrixIds.get(tf).forEach(matrixId -> add(() -> {
                    File logoFile = new File(tfDirectory, matrixId + ".svg");
                    IOUtils.copy(new URL(jasparLogoUrl.get() + matrixId + ".svg"), logoFile);

                    return true;
                }));
            });
        }};
    }
}
