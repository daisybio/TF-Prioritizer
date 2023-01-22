package org.exbio.tfprio.steps.logos;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import static org.exbio.pipejar.util.FileManagement.readLines;

public class FilterBindingSites extends ExecutableStep<Configs> {
    public final Map<String, OutputFile> outputFiles = new HashMap<>();
    private final Map<InputFile, OutputFile> bridge = new HashMap<>();

    public FilterBindingSites(Configs configs, Map<String, OutputFile> bindingSites) {
        super(configs, false, bindingSites.values());

        bindingSites.forEach((hm, hmDirectory) -> {
            OutputFile outputHm = addOutput(hm);
            outputFiles.put(hm, outputHm);
            bridge.put(addInput(hmDirectory), outputHm);
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((inputHm, outputHm) -> Arrays.stream(Objects.requireNonNull(inputHm.listFiles())).forEach(
                    tfInputFile -> add(() -> {
                        Map<Integer, List<String>> bindingSites =
                                readLines(tfInputFile).stream().collect(Collectors.groupingBy(String::length));

                        if (bindingSites.isEmpty()) {
                            return true;
                        }

                        List<String> mostBindingSites = bindingSites.entrySet().stream().max(
                                Comparator.comparingInt(map -> map.getValue().size())).get().getValue();

                        if (mostBindingSites.size() < 2) {
                            return true;
                        }

                        File tfOutputFile = new File(outputHm, tfInputFile.getName());

                        try (var writer = new BufferedWriter(new FileWriter(tfOutputFile))) {
                            mostBindingSites.forEach(line -> {
                                try {
                                    writer.write(line);
                                    writer.newLine();
                                } catch (IOException e) {
                                    throw new UncheckedIOException(e);
                                }
                            });
                        } catch (IOException e) {
                            throw new UncheckedIOException(e);
                        }

                        return true;
                    })));
        }};
    }
}
