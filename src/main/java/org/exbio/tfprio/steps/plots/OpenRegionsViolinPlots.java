package org.exbio.tfprio.steps.plots;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class OpenRegionsViolinPlots extends ExecutableStep<Configs> {
    public final OutputFile dataFile = addOutput("data.tsv");
    public final OutputFile plotFile = addOutput("plot.png");
    private final InputFile script;
    private final RequiredConfig<String> geneRegionFileName = new RequiredConfig<>(configs.tepic.geneRegionFileName);
    private final Map<String, Collection<InputFile>> hmDirectories = new HashMap<>();

    public OpenRegionsViolinPlots(Configs configs, Map<String, Map<String, Collection<OutputFile>>> tepicResults) {
        super(configs, false,
                tepicResults.values().stream().flatMap(x -> x.values().stream()).flatMap(Collection::stream).collect(
                        Collectors.toList()));

        tepicResults.forEach((group, hmMap) -> hmMap.forEach(
                (hm, files) -> hmDirectories.computeIfAbsent(hm, x -> new HashSet<>()).addAll(
                        files.stream().map(this::addInput).collect(Collectors.toSet()))));

        this.script =
                addInput(getClass().getResourceAsStream("openChromatinViolinPlots.R"), "openChromatinViolinPlots.R");
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                try (BufferedWriter writer = new BufferedWriter(new FileWriter(dataFile))) {
                    writer.write("HM\tLENGTH");
                    writer.newLine();

                    hmDirectories.entrySet().stream().map(dirEntry -> {
                        String hm = dirEntry.getKey();
                        Collection<File> files = dirEntry.getValue().stream().map(directory -> directory.listFiles(
                                file -> file.isFile() && file.getName().equals(geneRegionFileName.get()))).filter(
                                Objects::nonNull).map(foundFiles -> foundFiles[0]).collect(Collectors.toSet());
                        return new AbstractMap.SimpleEntry<>(hm, files);
                    }).forEach(hmFiles -> {
                        String hm = hmFiles.getKey();

                        hmFiles.getValue().forEach(file -> {
                            try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
                                reader.lines().skip(1).forEach(line -> {
                                    String regionString = line.split("\t")[0];
                                    String startEndString = regionString.substring(regionString.indexOf(":") + 1);
                                    String[] startEnd = startEndString.split("-");

                                    int start = Integer.parseInt(startEnd[0]);
                                    int end = Integer.parseInt(startEnd[1]);
                                    int length = end - start;

                                    try {
                                        writer.write(hm + "\t" + length);
                                        writer.newLine();
                                    } catch (IOException e) {
                                        throw new RuntimeException(e);
                                    }
                                });
                            } catch (IOException e) {
                                throw new RuntimeException(e);
                            }
                        });
                    });
                }

                String plotCommand = "Rscript " + script + " --input " + dataFile + " --output " + plotFile;

                executeAndWait(plotCommand, true);

                return true;
            });
        }};
    }
}
