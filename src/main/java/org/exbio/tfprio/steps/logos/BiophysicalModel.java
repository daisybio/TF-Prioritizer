package org.exbio.tfprio.steps.logos;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.exbio.pipejar.util.FileManagement.readLines;

public class BiophysicalModel extends ExecutableStep<Configs> {
    public final OutputFile outputFile = addOutput("models");
    private final InputFile dcgFile;
    private final InputFile pwmFile;
    private final RequiredConfig<File> pwmFileLocation = new RequiredConfig<>(configs.tepic.PWMs);

    public BiophysicalModel(Configs configs, OutputFile dcgFile) {
        super(configs, false, dcgFile);

        this.dcgFile = addInput(dcgFile);
        pwmFile = addInput(pwmFileLocation);
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

            Map<String, Collection<String>> tfBindingProfile = new HashMap<>();
            String currentTf = null;

            List<String> pwmLines;
            try {
                pwmLines = readLines(pwmFile);
            } catch (IOException e) {
                throw new RuntimeException(e);
            }

            for (String line : pwmLines) {
                if (line.startsWith("#")) {
                    continue;
                }

                if (line.startsWith(">")) {
                    currentTf = line.split("\t")[1];
                    continue;
                }

                if (tfNames.contains(currentTf)) {
                    tfBindingProfile.computeIfAbsent(currentTf, k -> new HashSet<>()).add(line);
                }
            }

            tfNames.forEach(tf -> add(() -> {
                List<String> bindingEnergies = Stream.of(tfBindingProfile.getOrDefault(tf, new HashSet<>()),
                        tf.contains("::") ? Arrays.stream(tf.split("::")).filter(tfBindingProfile::containsKey).flatMap(
                                singleTf -> tfBindingProfile.get(singleTf).stream()).toList() :
                                new HashSet<String>()).flatMap(Collection::stream).toList();

                File outputFile = new File(BiophysicalModel.this.outputFile, tf + ".tsv");

                try (var writer = new BufferedWriter(new FileWriter(outputFile))) {
                    bindingEnergies.forEach(line -> {
                        try {
                            writer.write(line);
                            writer.newLine();
                        } catch (IOException e) {
                            throw new RuntimeException(e);
                        }
                    });
                }


                return true;
            }));
        }};
    }
}
