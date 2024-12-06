package org.src;

import augmentedTree.IntervalTree;
import net.sf.samtools.AlignmentBlock;

import java.util.*;

public class Transcript {
    private final String transcriptId;
    private final int start;
    private final int stop;
    private final ArrayList<Exon> exonList = new ArrayList<>();
    private final char strand;
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
    }

    public String getTranscriptId() {
        return transcriptId;
    }

    public ArrayList<Exon> getExonList() {
        return this.exonList;
    }

    public ArrayList<Region> cut(int x1, int x2) {
        ArrayList<Region> cutRegions = new ArrayList<>();
        for (int i = 0; i < exonList.size(); i++) {
            Exon exon;
            if (strand == '-') {
                exon = exonList.get(exonList.size() - 1 - i);
            } else {
                exon = exonList.get(i);
            }

            // #----------------#
            //     <x1>----<x2>  → completely contained → add x1-x2 to set
            if (x1 >= exon.getStart() && x1 <= exon.getStop() && x2 >= exon.getStart() && x2 <= exon.getStop()) {
                cutRegions.add(new Region(x1, x2));
                break;
            }
            // #----------------#
            //     <x1>---------#
            else if (x1 >= exon.getStart() && x1 <= exon.getStop()) {
               cutRegions.add(new Region(x1, exon.getStop()));
            }
            // #----------------#
            // #----------------#
            else if (x1 <= exon.getStart() && x2 >= exon.getStop()) {
                cutRegions.add(new Region(exon.getStart(), exon.getStop()));
            }
            // #----------------#
            // #---------<x2>
            else if (x2 >= exon.getStart() && x2 <= exon.getStop()) {
                cutRegions.add(new Region(exon.getStart(), x2));
                break;
            }
        }
        return cutRegions;
    }

    public TreeSet<Region> cutSet(int x1, int x2) {
        // same as cut() but returns treeSet → useful for transcripotmic annotation
        TreeSet<Region> cutRegions = new TreeSet<>(
                Comparator.comparingInt(Region::getStart)
                        .thenComparingInt(Region::getStop)
        );
        for (int i = 0; i < exonList.size(); i++) {
            Exon exon;
            if (strand == '-') {
                exon = exonList.get(exonList.size() - 1 - i);
            } else {
                exon = exonList.get(i);
            }

            // #----------------#
            //     <x1>----<x2>  → completely contained → add x1-x2 to set
            if (x1 >= exon.getStart() && x1 <= exon.getStop() && x2 >= exon.getStart() && x2 <= exon.getStop()) {
                cutRegions.add(new Region(x1, x2));
                break;
            }
            // #----------------#
            //     <x1>---------#
            else if (x1 >= exon.getStart() && x1 <= exon.getStop()) {
                cutRegions.add(new Region(x1, exon.getStop()));
            }
            // #----------------#
            // #----------------#
            else if (x1 <= exon.getStart() && x2 >= exon.getStop()) {
                cutRegions.add(new Region(exon.getStart(), exon.getStop()));
            }
            // #----------------#
            // #---------<x2>
            else if (x2 >= exon.getStart() && x2 <= exon.getStop()) {
                cutRegions.add(new Region(exon.getStart(), x2));
                break;
            }
        }
        return cutRegions;
    }
}
