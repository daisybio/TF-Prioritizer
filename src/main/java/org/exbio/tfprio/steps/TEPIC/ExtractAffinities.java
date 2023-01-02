package org.exbio.tfprio.steps.TEPIC;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;

import java.io.File;
import java.util.*;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.FileManagement.makeSureFileExists;
import static org.exbio.pipejar.util.FileManagement.softLink;

public class ExtractAffinities extends ExecutableStep {
    public final Map<String, Map<String, Collection<OutputFile>>> outputFiles = new HashMap<>();
    private final Map<InputFile, OutputFile> bridge = new HashMap<>();

    public ExtractAffinities(Map<String, Map<String, Collection<OutputFile>>> tepicFiles) {
        super(false, tepicFiles.values().stream().flatMap(
                stringCollectionMap -> stringCollectionMap.values().stream().flatMap(Collection::stream)).toList());

        tepicFiles.forEach((group, hmMap) -> {
            OutputFile inputGroup = new OutputFile(inputDirectory, group);
            OutputFile outputGroup = new OutputFile(outputDirectory, group);

            hmMap.forEach((hm, sampleDirs) -> {
                OutputFile inputHm = new OutputFile(inputGroup, hm);
                OutputFile outputHm = new OutputFile(outputGroup, hm);

                sampleDirs.forEach(sampleDir -> {
                    InputFile inputDir = addInput(inputHm, sampleDir);
                    OutputFile outputFile = addOutput(outputHm, sampleDir.getName() + ".tsv");

                    bridge.put(inputDir, outputFile);
                    outputFiles.computeIfAbsent(group, s -> new HashMap<>()).computeIfAbsent(hm,
                            s -> new HashSet<>()).add(outputFile);
                });
            });
        });
    }

    @Override
    protected boolean doCreateFiles() {
        return false;
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((input, output) -> {
                add(() -> {
                    File[] allFiles = input.listFiles();

                    if (allFiles == null || allFiles.length == 0) {
                        return false;
                    }

                    List<String> attempts =
                            List.of("_Affinity_Gene_View_Filtered_TPM.txt", "_Affinity_Gene_View_Filtered.txt");

                    for (String attempt : attempts) {
                        List<File> matching =
                                Arrays.stream(allFiles).filter(file -> file.getName().endsWith(attempt)).toList();

                        if (matching.isEmpty()) {
                            continue;
                        }

                        if (matching.size() > 1) {
                            throw new IllegalStateException(
                                    "Found multiple files matching " + attempt + " in " + input);
                        }

                        File affinityFile = matching.get(0);

                        makeSureFileExists(output);
                        softLink(output, affinityFile);

                        return true;
                    }
                    logger.warn("Could not find affinity file for {}", input);
                    return false;
                });
            });
        }};
    }
}
