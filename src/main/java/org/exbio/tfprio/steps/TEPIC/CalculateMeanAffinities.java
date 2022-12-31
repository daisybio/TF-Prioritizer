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

public class CalculateMeanAffinities extends ExecutableStep {
    public final Map<String, Map<String, OutputFile>> outputFiles = new HashMap<>();
    private final Map<OutputFile, Collection<InputFile>> bridge = new HashMap<>();

    public CalculateMeanAffinities(Map<String, Map<String, Collection<OutputFile>>> tepicFiles) {
        super(false, tepicFiles.values().stream().flatMap(
                stringCollectionMap -> stringCollectionMap.values().stream()).flatMap(Collection::stream).toList());

        tepicFiles.forEach((group, hmMap) -> {
            OutputFile inGroup = new OutputFile(inputDirectory, group);
            OutputFile outGroup = new OutputFile(outputDirectory, group);
            hmMap.forEach((hm, sampleDirectories) -> {
                OutputFile outputFile = addOutput(outGroup, group + "_" + hm + ".tsv");

                outputFiles.computeIfAbsent(group, s -> new HashMap<>()).put(hm, outputFile);

                bridge.put(outputFile,
                        sampleDirectories.stream().map(sampleDirectory -> addInput(inGroup, sampleDirectory)).toList());
            });
        });
    }

    static void writeFile(Map<String, Map<String, Double>> geneTfAffinity, File outFile) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outFile))) {
            List<String> tfOrder = geneTfAffinity.values().stream().flatMap(
                    tfAffinity -> tfAffinity.keySet().stream()).distinct().sorted().toList();

            writer.write("Gene\t" + String.join("\t", tfOrder));
            writer.newLine();
            geneTfAffinity.forEach((gene, tfAffinity) -> {
                List<Double> affinityOrder =
                        tfOrder.stream().map(tf -> tfAffinity.getOrDefault(tf, Double.NaN)).toList();

                if (affinityOrder.stream().anyMatch(affinity -> affinity.isNaN())) {
                    return;
                }

                try {
                    writer.write(gene + "\t" +
                            affinityOrder.stream().map(Object::toString).collect(Collectors.joining("\t")));
                    writer.newLine();
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            });
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            // This step is not (yet) parallelized because it is not clear how to do it in a memory-efficient way
            add(() -> {
                bridge.forEach((outFile, sampleDirs) -> {
                    Map<String, Map<String, Double>> meanAffinity = getMeanAffinity(sampleDirs);
                    writeFile(meanAffinity, outFile);
                });

                return true;
            });
        }};
    }

    private Map<String, Map<String, Double>> getMeanAffinity(Collection<InputFile> sampleDirs) {
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
}
