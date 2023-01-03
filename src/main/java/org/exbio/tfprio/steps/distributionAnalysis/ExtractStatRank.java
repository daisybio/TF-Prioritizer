package org.exbio.tfprio.steps.distributionAnalysis;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.Callable;

public class ExtractStatRank extends ExecutableStep {
    public final Map<String, OutputFile> outputFiles = new HashMap<>();
    private final Map<InputFile, OutputFile> bridge = new HashMap<>();


    public ExtractStatRank(Map<String, OutputFile> hmStats) {
        super(false, hmStats.values());

        hmStats.forEach((hm, statsFile) -> {
            OutputFile output = addOutput(statsFile.getName());
            bridge.put(addInput(statsFile), output);
            outputFiles.put(hm, output);
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((input, output) -> add(() -> {
                try (BufferedReader reader = new BufferedReader(new FileReader(input));
                     BufferedWriter writer = new BufferedWriter(new FileWriter(output))) {
                    reader.readLine();
                    writer.write("TF\tRank");
                    writer.newLine();

                    reader.lines().forEachOrdered(line -> {
                        String[] split = line.split("\t");

                        if (split[1].equals("background")) {
                            return;
                        }

                        try {
                            writer.write(split[1] + "\t" + split[0]);
                            writer.newLine();
                        } catch (Exception e) {
                            e.printStackTrace();
                        }
                    });
                }

                return true;
            }));
        }};
    }
}
