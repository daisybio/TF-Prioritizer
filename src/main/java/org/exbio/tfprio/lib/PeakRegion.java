package org.exbio.tfprio.lib;

/**
 * Stores the score of a peak as well as its corresponding region
 */
public class PeakRegion extends Region {
    /**
     * The peak score
     */
    private final int score;
    /**
     * The line which was used in order to build this PeakRegion. Only set if the according constructor is used.
     */

    /**
     * Create a new PeakRegion based on its attributes.
     *
     * @param chromosome the chromosome the PeakRegion is located on
     * @param start      the start position of the region
     * @param end        the end position of the region
     * @param score      the peak score of this peak
     */
    public PeakRegion(String chromosome, int start, int end, int score) {
        super(chromosome, start, end);
        this.score = score;
    }


    /**
     * @return get the score of this PeakRegion
     */
    public int getScore() {
        return score;
    }

    public String toString() {
        return getChromosome() + "\t" + getStart() + "\t" + getEnd() + "\t" + getScore();
    }
}
