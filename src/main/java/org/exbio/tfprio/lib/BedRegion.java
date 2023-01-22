package org.exbio.tfprio.lib;

public class BedRegion extends Region {
    private final String signal;

    /**
     * Creates a new Region based on the given input data
     *
     * @param chromosome the chromosome this region is based on
     * @param start      the start position of this region
     * @param end        the end position of this region
     */
    public BedRegion(String chromosome, int start, int end, String signal) {
        super(chromosome, start, end);
        this.signal = signal;
    }

    private BedRegion(String[] split) {
        this(split[0], Integer.parseInt(split[1]), Integer.parseInt(split[2]), split[3]);
    }

    public BedRegion(String line) {
        this(line.split("\t"));
    }

    @Override
    public String toString() {
        return getChromosome() + "\t" + getStart() + "\t" + getEnd() + "\t" + signal;
    }
}
