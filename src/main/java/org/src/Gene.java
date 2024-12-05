package org.src;

import augmentedTree.Interval;
import augmentedTree.IntervalTree;
import net.sf.samtools.AlignmentBlock;

import java.util.*;

public class Gene implements Interval {
    private final int start;
    private final int end;
    private final String geneId;
    private final ArrayList<Transcript> transcriptList;
    private final HashMap<String, Transcript> transcriptMap;
    private final IntervalTree<Exon> exonTree;
    private final String geneName;
    private final String bioType;
    private final String chr;
    private final char strand;
    private ArrayList<Exon> combinedExonList;
    private TreeSet<Region> meltedRegions;

    private String sequence;


    public Gene(String geneId, int start, int end, String geneName, String chr, char strand, String bioType) {
        this.geneId = geneId;
        this.geneName = geneName;
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.strand = strand;
        this.bioType = bioType;
        this.transcriptList = new ArrayList<>();
        this.transcriptMap = new HashMap<>();
        this.exonTree = new IntervalTree<>();
    }

    public String getGeneId() {
        return geneId;
    }

    public void addTranscript(Transcript transcript) {
        transcriptList.add(transcript);
        transcriptMap.put(transcript.getTranscriptId(), transcript);
    }

    public ArrayList<Transcript> getTranscriptList() {
        return transcriptList;
    }

    public HashMap<String, Transcript> getTranscriptMap() {
        return transcriptMap;
    }

    public Transcript getLastTranscript() {
        if (!transcriptList.isEmpty()) {
            return transcriptList.get(transcriptList.size() - 1);
        }
        return null;
    }

    public String getChr() {
        return chr;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
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

    public char getStrand() {
        return strand;
    }

    public String getSeq() {
        return sequence;
    }

    public void addExon(int start, int end, int pos, String transcriptId) {
        this.exonTree.add(new Exon(start, end, pos, transcriptId));
    }

    public IntervalTree<Exon> getExonTree() {
        return exonTree;
    }

    public String getBioType() {
        return bioType;
    }

    public void initCombinedExonsList() {
        this.combinedExonList = new ArrayList<>();
        for (Transcript transcript : transcriptList) {
            for (int i = 0; i < transcript.getExonList().size(); i++) {
                Exon exon;
                if (strand == '-') {
                    exon = transcript.getExonList().get(transcript.getExonList().size() - 1 - i);
                } else {
                    exon = transcript.getExonList().get(i);
                }
                combinedExonList.add(exon);
            }
        }
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
}

