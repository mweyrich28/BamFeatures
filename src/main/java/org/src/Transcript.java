package org.src;

import augmentedTree.IntervalTree;

import java.util.*;

public class Transcript {
    private final String transcriptId;
    private final int start;
    private final int stop;
    private final ArrayList<Exon> exonList = new ArrayList<>();
    private final char strand;
    private final IntervalTree<Exon> exonTree = new IntervalTree<Exon>();

    private String transcriptSeq; // patched together using its exons

    public Transcript(String transcriptId, char strand, int start, int stop) {
        this.transcriptId = transcriptId;
        this.strand = strand;
        this.start = start;
        this.stop = stop;
    }

    public void addExon(int start, int end, int pos) {
        Exon exon = new Exon(start, end, pos, end - start + 1);
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

    public HashSet<String> cut(int x1, int x2) {
        HashSet<String> cutRegions = new HashSet<>();
        for (int i = 0; i < exonList.size(); i++) {
            Exon exon;
            if (strand == '-') {
                exon = exonList.get(exonList.size() - 1 - i);
            } else {
                exon = exonList.get(i);
            }

            // #----------------#
            //     x1----x2  → completely contained → add x1-x2 to set
            if (x1 >= exon.getStart() && x1 <= exon.getStop() && x2 >= exon.getStart() && x2 <= exon.getStop()) {
                cutRegions.add(x1 + "-" + x2);
                break;
            }
            else if (x1 >= exon.getStart() && x1 <= exon.getStop()) {
               cutRegions.add(x1 + "-" + exon.getStop());
            }
            else if (x1 <= exon.getStart() && x2 >= exon.getStop()) {
                cutRegions.add((exon.getStart() + 1) + "-" + (exon.getStop() - 1));
            }
            else if (x2 < exon.getStop()) {
                cutRegions.add(exon.getStart() + "-" + x2);
                break;
            }
        }
        return cutRegions;
    }
}
