package org.exbio.tfprio.steps.ehmm;

import org.apache.commons.lang3.builder.HashCodeBuilder;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.nio.file.Files;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class LearnBackgroundModel extends ExecutableStep<Configs> {
    private final OutputFile outputDir = addOutput("out");
    public final OutputFile outputFile = addOutput(outputDir,"BackgroundModel.RData");
    private final InputFile bedFile;
    private final InputFile bamDir;

    private final RequiredConfig<Integer> nStates = new RequiredConfig<>(configs.ehmm.nStates);
    private final RequiredConfig<Double> pseudoCount = new RequiredConfig<>(configs.ehmm.pseudoCount);
    private final Integer nBins = new RequiredConfig<>(configs.ehmm.nBins).get()*10;
    private final InputFile learnModelRscript;

    public LearnBackgroundModel(Configs configs, OutputFile bedFile, OutputFile bamDir) {
        super(configs, false, bamDir, bedFile);
        this.bedFile = addInput(bedFile);
        this.bamDir = addInput(bamDir);
        this.learnModelRscript = addInput(getClass().getResourceAsStream("learnModel.R"), "learnModel.R");
    }

    private static class Region implements Comparable<Region>{
        String chr;
        long start;
        long stop;
        char strand;
        String label;
        String score;

        public Region(String bedLine){
            String[] split = bedLine.split("\t");
            this.chr = split[0];
            this.start = Long.parseLong(split[1]);
            this.stop = Long.parseLong(split[2]);
            this.label = split[3];
            this.score = split.length < 5 ? "666": split[4];
            this.strand = split.length < 6 ? '*': split[5].charAt(0);
        }

        @Override
        public boolean equals(Object obj) {
            if (this.getClass() != obj.getClass()) return false;
            if (this == obj) return true;
            Region other = (Region) obj;
            return this.chr.equals(other.chr) &&
                    this.start == other.start &&
                    this.stop == other.stop;
        }

        @Override
        public int hashCode() {
            return new HashCodeBuilder()
                    .append(this.chr)
                    .append(this.start)
                    .append(this.stop)
                    .toHashCode();
        }

        @Override
        public int compareTo(Region o) {
            return this.chr.compareTo(o.chr)*10+(int) (this.start - o.start);
        }

        @Override
        public String toString(){
            return String.join("\t",
                    this.chr, String.valueOf(this.start), String.valueOf(this.stop),
                    this.label, this.score, String.valueOf(this.strand));
        }
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                // build background model
                String ehmmCommand = String.join(" ","Rscript",
                        learnModelRscript.getAbsolutePath(),
                        "-r", bedFile.getAbsolutePath(),
                        "-m", bamDir.getAbsolutePath(),
                        "-f", "BackgroundModel",
                        "-n", nStates.toString(),
                        "-b", nBins.toString(),
                        "-p", pseudoCount.toString(),
                        "-o", outputDir.getAbsolutePath());
                executeAndWait(ehmmCommand, true);
                return true;
            });
        }};
    }
}
