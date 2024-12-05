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

    @Override
    public String toString() {
        return getStart() + "-" + getStop();
    }

    @Override
    public int hashCode() {
        int result = 17;
        result = 31 * result + start;
        result = 31 * result + stop;
        return result;
    }
    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;
        Region region = (Region) obj;
        return start == region.start && stop == region.stop;
    }
    public boolean contains(Region other) {
        return this.start <= other.start && this.stop >= other.stop;
    }
}
