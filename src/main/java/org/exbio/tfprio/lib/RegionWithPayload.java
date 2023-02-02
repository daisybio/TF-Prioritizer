package org.exbio.tfprio.lib;

public class RegionWithPayload<T> extends Region {
    private final T payload;

    public RegionWithPayload(String chromosome, int start, int end, T payload) {
        super(chromosome, start, end);
        this.payload = payload;
    }

    public T getPayload() {
        return payload;
    }
}
