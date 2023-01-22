package org.exbio.tfprio.lib;

import java.util.HashSet;
import java.util.Set;

public class MultiGeneRegion extends Region {
    private final Set<String> ids;

    public MultiGeneRegion(String chromosome, int start, int end, Set<String> ids) {
        super(chromosome, start, end);
        this.ids = ids;
    }

    public MultiGeneRegion(String chromosome, int start, int end, String id) {
        super(chromosome, start, end);
        this.ids = new HashSet<>() {{
            add(id);
        }};
    }

    public void addId(String id) {
        ids.add(id);
    }

    public Set<String> getIds() {
        return ids;
    }

    public Region getRegion() {
        return new Region(getChromosome(), getStart(), getEnd());
    }

    public String toString() {
        return getChromosome() + "\t" + getStart() + "\t" + getEnd() + "\t" + (String.join(";", getIds()));
    }
}
