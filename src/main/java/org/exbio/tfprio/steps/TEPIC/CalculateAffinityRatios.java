package org.exbio.tfprio.steps.TEPIC;

import org.apache.commons.lang3.tuple.Pair;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.Callable;

import static org.exbio.tfprio.steps.TEPIC.CalculateMeanAffinities.writeFile;

public class CalculateAffinityRatios extends ExecutableStep<Configs> {
    public final Map<String, Map<String, OutputFile>> outputFiles = new HashMap<>();
    private final Map<OutputFile, Pair<InputFile, InputFile>> bridge = new HashMap<>();

    public CalculateAffinityRatios(Configs configs, Map<String, Map<String, OutputFile>> meanAffinities) {
        super(configs, false, meanAffinities.values().stream().flatMap(
                stringCollectionMap -> stringCollectionMap.values().stream()).toList());

        meanAffinities.forEach((group1, hmMap1) -> meanAffinities.forEach((group2, hmMap2) -> {
            if (group1.compareTo(group2) >= 0) {
                return;
            }

            String pairing = group1 + "_" + group2;

            OutputFile outPairing = new OutputFile(outputDirectory, pairing);
            OutputFile inPairing = new OutputFile(inputDirectory, pairing);

            outputFiles.computeIfAbsent(pairing, s -> new HashMap<>());

            hmMap1.forEach((hm, file1) -> {
                if (!hmMap2.containsKey(hm)) {
                    return;
                }
                OutputFile file2 = hmMap2.get(hm);
                OutputFile outHm = addOutput(outPairing, hm + ".tsv");

                bridge.put(outHm, Pair.of(addInput(inPairing, file1), addInput(inPairing, file2)));
                outputFiles.get(pairing).put(hm, outHm);
            });
        }));
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            // This step is not (yet) parallelized because it is not clear how to do it in a memory-efficient way
            add(() -> {
                bridge.forEach((outputFile, pair) -> {
                    Map<String, Map<String, Double>> meanAffinity1 = readAffinityFile(pair.getLeft());
                    Map<String, Map<String, Double>> meanAffinity2 = readAffinityFile(pair.getRight());

                    Map<String, Map<String, Double>> affinityRatios =
                            meanAffinity1.entrySet().stream().map(geneEntry -> {
                                String gene = geneEntry.getKey();
                                Map<String, Double> tfAffinity1 = geneEntry.getValue();
                                Map<String, Double> tfAffinity2 = meanAffinity2.getOrDefault(gene, new HashMap<>(0));

                                Map<String, Double> tfAffinityRatios = tfAffinity1.entrySet().stream().map(tfEntry -> {
                                    String tf = tfEntry.getKey();
                                    double affinity1 = tfEntry.getValue();

                                    if (!tfAffinity2.containsKey(tf)) {
                                        return Pair.of(tf, Double.NaN);
                                    }

                                    double affinity2 = tfAffinity2.get(tf);

                                    if (affinity2 == 0) {
                                        if (affinity1 == 0) {
                                            return Pair.of(tf, 1.0);
                                        } else {
                                            return Pair.of(tf, Double.NaN);
                                        }
                                    }

                                    return Pair.of(tf, affinity1 / affinity2);
                                }).collect(HashMap::new, (map, p) -> map.put(p.getLeft(), p.getRight()),
                                        HashMap::putAll);
                                return Pair.of(gene, tfAffinityRatios);
                            }).collect(HashMap::new, (map, p) -> map.put(p.getLeft(), p.getRight()), HashMap::putAll);

                    writeFile(affinityRatios, outputFile);
                });
                return true;
            });
        }};
    }

    private Map<String, Map<String, Double>> readAffinityFile(File file) {
        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String[] tfLocations = br.readLine().split("\t");

            return br.lines().map(line -> line.split("\t")).map(line -> {
                Map<String, Double> tfAffinity = new HashMap<>();
                String gene = line[0];
                for (int i = 1; i < line.length; i++) {
                    String tf = tfLocations[i];
                    double affinity = Double.parseDouble(line[i]);

                    tfAffinity.put(tf, affinity);
                }
                return Pair.of(gene, tfAffinity);
            }).collect(HashMap::new, (map, p) -> map.put(p.getLeft(), p.getRight()), HashMap::putAll);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
