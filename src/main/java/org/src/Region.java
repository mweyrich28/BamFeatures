package org.src;

import augmentedTree.Interval;

public class Region implements Interval {
    private int start;
    private int stop;

    public Region(int start, int end) {
        this.start = start;
        this.stop = end;
    }

    @Override
    public int getStart() {
        return start;
    }

    @Override
    public int getStop() {
        return stop;
    }

    public void setStop(int stop) {
        this.stop = stop;
    }

    public void setStart(int start) {
        this.start = start;
    }
}
