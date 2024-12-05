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

    public IntervalTree<Exon> getExonTree() {
        return exonTree;
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
            //     x1----x2  → completely contained → add x1-x2 to set
            if (x1 >= exon.getStart() && x1 <= exon.getStop() && x2 >= exon.getStart() && x2 <= exon.getStop()) {
                cutRegions.add(new Region(x1, x2));
                break;
            }
            else if (x1 >= exon.getStart() && x1 <= exon.getStop()) {
               cutRegions.add(new Region(x1, exon.getStop()));
            }
            else if (x1 <= exon.getStart() && x2 >= exon.getStop()) {
                cutRegions.add(new Region(exon.getStart(), exon.getStop()));
            }
            else if (x2 < exon.getStop()) {
                cutRegions.add(new Region(exon.getStart(), x2));
                break;
            }
        }
        return cutRegions;
    }

    public HashSet<Region> meltRegions(ArrayList<Region> regions) {
        HashSet<Region> melted = new HashSet<>();

        if (regions.size() == 1) {
            melted.addAll(regions);
            return melted;
        } else {
            // Sort regions by start position
            Collections.sort(regions, Comparator.comparingInt(Region::getStart));

            // Initialize the first region to merge
            Region current = new Region(regions.get(0).getStart(), regions.get(0).getStop());

            for (int i = 1; i < regions.size(); i++) {
                Region exon = regions.get(i);

                // Merge overlapping or adjacent regions
                if (exon.getStart() <= current.getStop() + 1) {
                    current.setStop(Math.max(current.getStop(), exon.getStop()));
                } else {
                    // Add the merged region to the set
                    melted.add(current);

                    // Start a new region to merge
                    current = new Region(exon.getStart(), exon.getStop());
                }
            }

            // Add the last region
            melted.add(current);
            return melted;
        }
    }
}
