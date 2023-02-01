package org.exbio.tfprio.steps.rnaSeq;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.FileManagement.readLines;

public class FilterENdb extends ExecutableStep<Configs> {
    public final OutputFile outputFile;
    private final InputFile enDBfile;
    private final RequiredConfig<String> species = new RequiredConfig<>(configs.inputConfigs.genome);

    public FilterENdb(Configs configs) {
        super(configs);
        enDBfile = addInput(getClass().getResourceAsStream("ENdb.bed"), "ENdb.bed");
        outputFile = addOutput("out.bed");
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                String speciesName = species.get();
                String symbol = speciesName.replaceAll("\\d", "");

                try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
                    readLines(enDBfile).stream().skip(1).map(line -> line.split("\t")).filter(
                            splittedLine -> splittedLine[4].replaceAll("\\d", "").equals(symbol)).forEach(line -> {
                        try {
                            writer.write(String.join("\t", line));
                            writer.newLine();
                        } catch (IOException e) {
                            throw new RuntimeException(e);
                        }
                    });
                }
                return true;
            });
        }};
    }
}
