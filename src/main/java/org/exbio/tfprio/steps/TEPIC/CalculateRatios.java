package org.exbio.tfprio.steps.TEPIC;

import org.apache.commons.lang3.tuple.Pair;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class CalculateRatios extends ExecutableStep {
    public final Collection<OutputFile> outputFiles = new HashSet<>();
    private final Map<OutputFile, Pair<Collection<InputFile>, Collection<InputFile>>> bridge = new HashMap<>();

    public CalculateRatios(Map<String, Map<String, Collection<OutputFile>>> tepicFiles) {
        super(false, tepicFiles.values().stream().flatMap(
                stringCollectionMap -> stringCollectionMap.values().stream()).flatMap(Collection::stream).toList());

        tepicFiles.forEach((group1, hmMap1) -> tepicFiles.forEach((group2, hmMap2) -> {
            if (group1.compareTo(group2) >= 0) {
                return;
            }

            String pairing = group1 + "_" + group2;

            OutputFile outPairing = new OutputFile(outputDirectory, pairing);
            OutputFile inPairing = new OutputFile(inputDirectory, pairing);

            hmMap1.forEach((hm, files1) -> {
                if (!hmMap2.containsKey(hm)) {
                    return;
                }

                Collection<OutputFile> files2 = hmMap2.get(hm);

                OutputFile outHm = addOutput(outPairing, hm);
                OutputFile inHm = new OutputFile(inPairing, hm);

                OutputFile group1Dir = new OutputFile(inHm, group1);
                OutputFile group2Dir = new OutputFile(inHm, group2);

                Collection<InputFile> sampleDirs1 =
                        files1.stream().map(outputFile -> addInput(group1Dir, outputFile)).toList();
                Collection<InputFile> sampleDirs2 =
                        files2.stream().map(outputFile -> addInput(group2Dir, outputFile)).toList();

                bridge.put(outHm, Pair.of(sampleDirs1, sampleDirs2));
            });
        }));
    }

    private static Map<String, Map<String, Double>> getMeanAffinity(Collection<InputFile> sampleDirs) {
        return sampleDirs.stream()
                         // Get affinity files
                         .map(sampleDir -> {
                             File[] consideredFiles = sampleDir.listFiles(
                                     file -> file.isFile() && file.getName().contains("Affinity_Gene_View_Filtered"));

                             if (consideredFiles == null || consideredFiles.length != 1) {
                                 throw new RuntimeException(
                                         "Expected exactly one file in " + sampleDir + " but found " +
                                                 (consideredFiles == null ? 0 : consideredFiles.length));
                             }

                             return consideredFiles[0];
                         })
                         // Get affinities
                         .map(affinityFile -> {
                             Map<String, Map<String, Double>> geneTfAffinity;

                             try (BufferedReader reader = new BufferedReader(new FileReader(affinityFile))) {
                                 String[] tfLocations = reader.readLine().split("\t");

                                 geneTfAffinity = reader.lines().map(line -> {
                                     String[] split = line.split("\t");
                                     String gene = split[0];
                                     Map<String, Double> tfAffinity = IntStream.range(1, split.length).mapToObj(
                                             index -> Pair.of(tfLocations[index],
                                                     Double.parseDouble(split[index]))).collect(HashMap::new,
                                             (map, p) -> map.put(p.getLeft(), p.getRight()), HashMap::putAll);
                                     return Pair.of(gene, tfAffinity);
                                 }).collect(HashMap::new, (map, p) -> map.put(p.getLeft(), p.getRight()),
                                         HashMap::putAll);
                             } catch (IOException e) {
                                 throw new RuntimeException(e);
                             }

                             return geneTfAffinity;
                         })
                         // Concatenate affinities
                         .reduce(new HashMap<String, Map<String, List<Double>>>(), (map, geneTfAffinity) -> {
                             geneTfAffinity.forEach((gene, tfAffinity) -> tfAffinity.forEach(
                                     (tf, affinity) -> map.computeIfAbsent(gene, s -> new HashMap<>()).computeIfAbsent(
                                             tf, s -> new ArrayList<>()).add(affinity)));
                             return map;
                         }, (map1, map2) -> {
                             map2.forEach((gene, tfAffinity) -> tfAffinity.forEach(
                                     (tf, newAffinities) -> map1.computeIfAbsent(gene,
                                             s -> new HashMap<>()).computeIfAbsent(tf, s -> new ArrayList<>()).addAll(
                                             newAffinities)));
                             return map1;
                         }).entrySet().stream()
                         // Average affinities
                         .map(entry -> {
                             Map<String, Double> tfAffinity = entry.getValue().entrySet().stream().map(
                                     tfEntry -> Pair.of(tfEntry.getKey(), tfEntry.getValue().stream().mapToDouble(
                                             Double::doubleValue).average().orElseThrow())).collect(HashMap::new,
                                     (map, p) -> map.put(p.getLeft(), p.getRight()), HashMap::putAll);
                             return Pair.of(entry.getKey(), tfAffinity);
                         }).collect(HashMap::new, (map, p) -> map.put(p.getLeft(), p.getRight()), HashMap::putAll);
    }

    @Override
    protected boolean mayBeSkipped() {
        return false;
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            // This step is not (yet) parallelized because it is not clear how to do it in a memory-efficient way
            add(() -> {
                bridge.forEach((outDir, pair) -> {
                    Collection<InputFile> sampleDirs1 = pair.getLeft();
                    Collection<InputFile> sampleDirs2 = pair.getRight();

                    Map<String, Map<String, Double>> meanAffinity1 = getMeanAffinity(sampleDirs1);
                    writeFile(meanAffinity1, new File(outDir, "mean_affinity_1.tsv"));

                    Map<String, Map<String, Double>> meanAffinity2 = getMeanAffinity(sampleDirs2);
                    writeFile(meanAffinity2, new File(outDir, "mean_affinity_2.tsv"));

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
                                            return Pair.of(tf, Double.NaN);
                                        } else {
                                            return Pair.of(tf, Double.POSITIVE_INFINITY);
                                        }
                                    }

                                    return Pair.of(tf, affinity1 / affinity2);
                                }).collect(HashMap::new, (map, p) -> map.put(p.getLeft(), p.getRight()),
                                        HashMap::putAll);
                                return Pair.of(gene, tfAffinityRatios);
                            }).collect(HashMap::new, (map, p) -> map.put(p.getLeft(), p.getRight()), HashMap::putAll);

                    writeFile(affinityRatios, new File(outDir, "affinity_ratios.tsv"));
                });
                return true;
            });
        }};
    }

    private void writeFile(Map<String, Map<String, Double>> geneTfAffinity, File outFile) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outFile))) {
            writer.write("Gene\t" + geneTfAffinity.values().stream().flatMap(
                    tfAffinity -> tfAffinity.keySet().stream()).distinct().sorted().collect(Collectors.joining("\t")));
            writer.newLine();
            geneTfAffinity.forEach((gene, tfAffinity) -> {
                try {
                    writer.write(gene + "\t" +
                            tfAffinity.values().stream().map(Object::toString).collect(Collectors.joining("\t")));
                    writer.newLine();
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            });
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
