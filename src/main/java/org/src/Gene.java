package org.src;

import augmentedTree.Interval;
import augmentedTree.IntervalTree;
import net.sf.samtools.AlignmentBlock;

import java.util.*;

public class Gene implements Interval {
    private final int start;
    private final int end;
    private int meltedLength = 0;

    private final String geneId;
    private final ArrayList<Transcript> transcriptList;
    private final String geneName;
    private final String bioType;
    private final String chr;
    private final char strand;
    private TreeSet<Region> meltedRegions;
    public Gene(String geneId, int start, int end, String geneName, String chr, char strand, String bioType) {
        this.geneId = geneId;
        this.geneName = geneName;
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.strand = strand;
        this.bioType = bioType;
        this.transcriptList = new ArrayList<>();
    }

    public String getGeneId() {
        return geneId;
    }

    public void addTranscript(Transcript transcript) {
        transcriptList.add(transcript);
    }

    public ArrayList<Transcript> getTranscriptList() {
        return transcriptList;
    }

    public String getChr() {
        return chr;
    }

    @Override
    public int getStart() {
        return start;
    }

    @Override
    public int getStop() {
        return end;
    }

    public int getEnd() {
        return end;
    }
    public Transcript getLastTranscript() {
        if (!transcriptList.isEmpty()) {
            return transcriptList.get(transcriptList.size() - 1);
        }
        return null;
    }

    public char getStrand() {
        return strand;
    }

    public String getBioType() {
        return bioType;
    }

    public void melt() {

        ArrayList<Exon> allExons = new ArrayList<>();
        for (Transcript transcript: transcriptList) {
            allExons.addAll(transcript.getExonList());
        }

        Collections.sort(allExons, Comparator.comparingInt(Exon::getStart));
        TreeSet<Region> meltedRegions = new TreeSet<>(
                Comparator.comparingInt(Region::getStart)
                        .thenComparingInt(Region::getStop)
        );


        if (!allExons.isEmpty()) {
            Exon first = allExons.getFirst();
            Region current = new Region(first.getStart(), first.getStop());

            for (int i = 1; i < allExons.size(); i++) {
                Exon exon = allExons.get(i);

                if (exon.getStart() <= current.getStop() + 1) {
                    current.setStop(Math.max(current.getStop(), exon.getStop()));
                } else {
                    meltedRegions.add(current);
                    current = new Region(exon.getStart(), exon.getStop());
                }
            }

            meltedRegions.add(current);
        }
        this.meltedRegions = meltedRegions;
    }

    public TreeSet<Region> getMeltedRegions() {
        if (this.meltedRegions == null) {
            melt();
        }
        return this.meltedRegions;
    }

    public int getMeltedLength() {
        // make sure to only calculate once
        if (this.meltedLength == 0) {
            TreeSet<Region> meltedRegions = getMeltedRegions();
            int length = 0;
            for (Region region : meltedRegions) {
                length += region.getStop() - region.getStart() - 1;
            }
            this.meltedLength = length;
        }
        return meltedLength;
    }
}

