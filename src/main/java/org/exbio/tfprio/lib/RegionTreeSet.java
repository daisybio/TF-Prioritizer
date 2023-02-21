package org.exbio.tfprio.lib;

import java.util.*;

public class RegionTreeSet extends TreeSet<Region> {
    public RegionTreeSet(Collection<Region> regions) {
        super();
        addMergeAll(regions);
    }

    public boolean addMergeAll(Collection<Region> regions) {
        return regions.stream().map(this::addMerge).reduce(Boolean::logicalOr).orElse(false);
    }

    public boolean addMerge(Region region) {
        int start = region.getStart();
        int end = region.getEnd();
        String chromosome = region.getChromosome();

        Set<Region> overlapping = new HashSet<>() {{
            addAll(getOverlapping(region));
        }};

        if (overlapping.isEmpty()) {
            return add(region);
        } else {
            int newStart =
                    Math.min(start, overlapping.stream().map(Region::getStart).min(Integer::compareTo).orElseThrow());
            int newEnd = Math.max(end, overlapping.stream().map(Region::getEnd).max(Integer::compareTo).orElseThrow());
            Region addRegion = new Region(chromosome, newStart, newEnd);

            removeAll(overlapping);
            return add(addRegion);
        }
    }

    public SortedSet<Region> getOverlapping(Region region) {
        Region startRegion = new Region(region.getChromosome(), region.getStart(), region.getStart());
        Region endRegion = new Region(region.getChromosome(), region.getEnd(), region.getEnd());

        return subSet(startRegion, true, endRegion, true);
    }

    public boolean hasOverlap(Region region) {
        return !getOverlapping(region).isEmpty();
    }

    public boolean hasOverlap(Region region, int margin) {
        return hasOverlap(new Region(region.getChromosome(), region.getStart() - margin, region.getEnd() + margin));
    }
}
