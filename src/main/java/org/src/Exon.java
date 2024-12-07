package org.src;

public class Exon {
    private int length;
    private final int start;
    private final int stop;
    private int pos;
    private String transcriptId;

    public Exon(int start, int end, int pos, int length) {
        this.length = length;
        this.start = start;
        this.stop = end;
        this.pos = pos;
    }

    public Exon(int start, int end, int pos, String tId) {
        this.start = start;
        this.stop = end;
        this.pos = pos;
        this.length = stop - start + 1;
        this.transcriptId = tId;
    }

    public int getStart() {
        return start;
    }

    public int getStop() {
        return stop;
    }

    public int getPos() {
        return pos;
    }

    public void setPos(int pos) {
        this.pos = pos;
    }

    @Override
    public String toString() {
        return this.transcriptId + " " + this.start + "-" + this.stop + " " + "[" + this.pos +"] Length:" + this.length;
    }

    public int getLength() {
        return length;
    }
}
