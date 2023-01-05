package org.exbio.tfprio.steps.TEPIC;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;
import org.exbio.tfprio.lib.Region;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

public class CreateBindingRegionsBedFiles extends ExecutableStep {
    public final Map<String, Map<String, OutputFile>> outputFiles = new HashMap<>();
    private final Map<OutputFile, Collection<InputFile>> bridge = new HashMap<>();
    private final RequiredConfig<Double> affinityCutoff = new RequiredConfig<>(Configs.tepic.affinityCutoff);
    private final RequiredConfig<String> sequenceFileName = new RequiredConfig<>(Configs.tepic.sequenceFileName);

    private final Map<OutputFile, String> outputDirectoryGroup = new HashMap<>();

    public CreateBindingRegionsBedFiles(Map<String, Map<String, Collection<OutputFile>>> tepicResults) {
        super(false,
                tepicResults.values().stream().flatMap(x -> x.values().stream()).flatMap(Collection::stream).collect(
                        Collectors.toList()));
        tepicResults.forEach((group, hmMap) -> {
            outputFiles.put(group, new HashMap<>());
            OutputFile outputGroup = new OutputFile(outputDirectory, group);
            hmMap.forEach((hm, sampleDirectories) -> {
                OutputFile outputDirectory = addOutput(outputGroup, hm);
                outputFiles.get(group).put(hm, outputDirectory);
                outputDirectoryGroup.put(outputDirectory, group);
                bridge.put(outputDirectory, sampleDirectories.stream().map(this::addInput).collect(Collectors.toSet()));
            });
        });
    }

    private static void addRegion(TreeSet<Region> a, Region b) {
        Region upper = a.ceiling(b);
        Region lower = a.floor(b);

        if (upper != null && upper.overlaps(b)) {
            upper.merge(b);
        } else if (lower != null && lower.overlaps(b)) {
            lower.merge(b);
        } else {
            a.add(b);
        }
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((outputDirectory, inputDirectories) -> add(() -> {
                Map<String, TreeSet<Region>> allRegions = inputDirectories.stream().map(
                                                                                  inputDirectory -> new File(inputDirectory, sequenceFileName.get())).map(sequenceFile -> {
                                                                              Map<String, TreeSet<Region>> sampleRegions = new HashMap<>();

                                                                              // Read sequence file
                                                                              try (BufferedReader reader = new BufferedReader(new FileReader(sequenceFile))) {
                                                                                  for (String line = reader.readLine(); line != null; line = reader.readLine()) {
                                                                                      if (line.startsWith("#") || line.startsWith("TF\t")) {
                                                                                          continue;
                                                                                      }

                                                                                      String[] split = line.split("\t");

                                                                                      // Regex to remove version numbers (e.g. ATF3(MA0605.1))
                                                                                      String tf = split[0].replaceAll("\\(.*?\\)", "");

                                                                                      double affinity = Double.parseDouble(split[1]);
                                                                                      int start = Integer.parseInt(split[3]);
                                                                                      int length = Integer.parseInt(split[4]);

                                                                                      if (affinity < affinityCutoff.get()) {
                                                                                          continue;
                                                                                      }

                                                                                      String area = split[6];
                                                                                      String chromosome = area.substring(1, area.indexOf(':'));
                                                                                      int areaStart = Integer.parseInt(area.substring(area.indexOf(':') + 1, area.indexOf('-')));
                                                                                      int areaEnd = Integer.parseInt(area.substring(area.indexOf('-') + 1));

                                                                                      int predictedStart = areaStart + start;
                                                                                      int predictedEnd = predictedStart + length;

                                                                                      sampleRegions.computeIfAbsent(tf, x -> new TreeSet<>()).add(
                                                                                              new Region(chromosome, predictedStart, predictedEnd));
                                                                                  }
                                                                              } catch (IOException e) {
                                                                                  throw new RuntimeException(e);
                                                                              }
                                                                              return sampleRegions;
                                                                          })
                                                                          // Collect regions from all samples to single map
                                                                          .reduce(new HashMap<>(), (a, b) -> {
                                                                              b.forEach(
                                                                                      (tf, regions) -> a.computeIfAbsent(
                                                                                              tf,
                                                                                              x -> new TreeSet<>()).addAll(
                                                                                              regions));
                                                                              return a;
                                                                          });

                Map<String, Collection<Region>> filteredRegions = allRegions.entrySet().stream()
                                                                            // Merge overlapping regions
                                                                            .map(entry -> {
                                                                                String tf = entry.getKey();
                                                                                TreeSet<Region> regions =
                                                                                        entry.getValue();

                                                                                TreeSet<Region> mergedRegions =
                                                                                        regions.stream().reduce(
                                                                                                new TreeSet<>(),
                                                                                                (a, b) -> {
                                                                                                    if (a.isEmpty()) {
                                                                                                        a.add(b);
                                                                                                        return a;
                                                                                                    }

                                                                                                    addRegion(a, b);
                                                                                                    return a;
                                                                                                }, (a, b) -> {
                                                                                                    b.forEach(
                                                                                                            region -> addRegion(
                                                                                                                    a,
                                                                                                                    region));
                                                                                                    return a;
                                                                                                });

                                                                                return new AbstractMap.SimpleEntry<>(tf,
                                                                                        mergedRegions);
                                                                            }).collect(
                                Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));

                filteredRegions.forEach((tf, regions) -> {
                    File outputFile = new File(outputDirectory, tf + ".bed");

                    writeRegionsToFile(tf, outputDirectoryGroup.get(outputDirectory), regions, outputFile);

                    if (tf.contains("::")) {
                        String[] tfParts = tf.split("::");

                        TreeSet<Region> merged =
                                Arrays.stream(tfParts).map(filteredRegions::get).filter(Objects::nonNull).flatMap(
                                        Collection::stream).reduce(new TreeSet<>(), (a, b) -> {
                                    if (a.isEmpty()) {
                                        a.add(b);
                                        return a;
                                    }

                                    addRegion(a, b);
                                    return a;
                                }, (a, b) -> {
                                    b.forEach(region -> addRegion(a, region));
                                    return a;
                                });
                        writeRegionsToFile(tf, outputDirectoryGroup.get(outputDirectory), merged,
                                new File(outputDirectory, tf + "_merged.bed"));
                    }
                });

                return true;
            }));
        }};
    }

    private void writeRegionsToFile(String tf, String group, Collection<Region> regions, File outputFile) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
            writer.write("track name=\"" + tf + "(@PRED_" + group + ")\"");
            writer.newLine();

            for (Region region : regions) {
                writer.write(region.getChromosome() + "\t" + region.getStart() + "\t" + region.getEnd());
                writer.newLine();
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
