package org.src;

import augmentedTree.IntervalTree;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.TreeMap;

public class Transcript {
    private final String transcriptId;
    private final int start;
    private final int stop;
    private final ArrayList<Exon> exonList = new ArrayList<>();
    private final String transcriptType;
    private final IntervalTree<Exon> exonTree = new IntervalTree<Exon>();

    private String transcriptSeq; // patched together using its exons

    public Transcript(String transcriptId, String transcriptType, int start, int stop) {
        this.transcriptId = transcriptId;
        this.transcriptType = transcriptType;
        this.start = start;
        this.stop = stop;
    }

    public void addExon(int start, int end, int pos) {
        Exon exon = new Exon(start, end, pos, end-start + 1);
        exonList.add(exon);
        exonTree.add(exon);
    }

    public String getTranscriptId() {
        return transcriptId;
    }

    public ArrayList<Exon> getExonList() {
        return this.exonList;
    }

    public void reversCdsList() {
        Collections.reverse(this.exonList);
    }

    public void setTranscriptSeq(String transcriptSeq) {
        this.transcriptSeq = transcriptSeq;
    }

    public String getTranscriptSeq() {
        return this.transcriptSeq;
    }

    public IntervalTree<Exon> getExonTree() {
        return exonTree;
    }
}
