package org.exbio.tfprio.lib;

import java.util.Collection;
import java.util.stream.Collectors;

public class BroadPeakRegion extends Region {
    private final double score, signalValue, pValue, qValue;
    private final String name;
    private final String strand;

    public BroadPeakRegion(String chromosome, int start, int end, String name, double score, String strand,
                           double signalValue, double pValue, double qValue) {
        super(chromosome, start, end);
        this.score = score;
        this.signalValue = signalValue;
        this.pValue = pValue;
        this.qValue = qValue;
        this.name = name;
        this.strand = strand;
    }

    private BroadPeakRegion(String[] splittedLine) {
        this(splittedLine[0], Integer.parseInt(splittedLine[1]), Integer.parseInt(splittedLine[2]), splittedLine[3],
                Double.parseDouble(splittedLine[4]), splittedLine[5], Double.parseDouble(splittedLine[6]),
                Double.parseDouble(splittedLine[7]), Double.parseDouble(splittedLine[8]));
    }

    public BroadPeakRegion(String line) {
        this(line.split("\t"));
    }

    public static BroadPeakRegion merge(Collection<BroadPeakRegion> overlapping) {
        if (overlapping.isEmpty()) {
            throw new IllegalArgumentException("Cannot merge empty collection");
        }

        int start = overlapping.stream().map(Region::getStart).min(Integer::compareTo).get();
        int end = overlapping.stream().map(Region::getEnd).max(Integer::compareTo).get();
        double score = overlapping.stream().mapToDouble(BroadPeakRegion::getScore).average().getAsDouble();
        double signalValue = overlapping.stream().mapToDouble(BroadPeakRegion::getSignalValue).average().getAsDouble();
        double pValue = overlapping.stream().mapToDouble(BroadPeakRegion::getpValue).average().getAsDouble();
        double qValue = overlapping.stream().mapToDouble(BroadPeakRegion::getqValue).average().getAsDouble();
        String name = overlapping.stream().map(BroadPeakRegion::getName).collect(Collectors.joining(","));
        String strand = ".";
        return new BroadPeakRegion(overlapping.iterator().next().getChromosome(), start, end, name, score, strand,
                signalValue, pValue, qValue);
    }

    public static BroadPeakRegion between(BroadPeakRegion a, BroadPeakRegion b) {
        if (!a.getChromosome().equals(b.getChromosome())) {
            throw new IllegalArgumentException("Cannot merge regions on different chromosomes");
        }

        int start = a.getEnd();
        int end = b.getStart();
        double score = (a.getScore() + b.getScore()) / 2;
        double signalValue = (a.getSignalValue() + b.getSignalValue()) / 2;
        double pValue = (a.getpValue() + b.getpValue()) / 2;
        double qValue = (a.getqValue() + b.getqValue()) / 2;
        String name = a.getName() + "," + b.getName();
        String strand = ".";
        return new BroadPeakRegion(a.getChromosome(), start, end, name, score, strand, signalValue, pValue, qValue);
    }

    @Override
    public String toString() {
        return getChromosome() + "\t" + getStart() + "\t" + getEnd() + "\t" + name + "\t" + score + "\t" + strand +
                "\t" + signalValue + "\t" + pValue + "\t" + qValue;
    }

    public double getScore() {
        return score;
    }

    public double getSignalValue() {
        return signalValue;
    }

    public double getpValue() {
        return pValue;
    }

    public double getqValue() {
        return qValue;
    }

    public String getName() {
        return name;
    }

    public String getStrand() {
        return strand;
    }
}
