package org.src;

import augmentedTree.Interval;

public class Intron implements Interval {
    private final int start;
    private final int stop;
    private final boolean fw;

    public Intron(int start, int stop, boolean fw) {
        this.start = start;
        this.stop = stop;
        this.fw = fw;
    }

    @Override
    public int getStart() {
        return start;
    }

    @Override
    public int getStop() {
        return stop;
    }

    public boolean isFw() {
        return fw;
    }
}

