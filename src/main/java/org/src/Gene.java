package org.src;

import augmentedTree.Interval;
import augmentedTree.IntervalTree;
import net.sf.samtools.AlignmentBlock;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

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

    public void invertTranscripts() {
        for (int i = 0; i < transcriptList.size(); i++) {
            Transcript currTranscript = transcriptList.get(i);
            currTranscript.reversCdsList();
            for (int j = 0; j < currTranscript.getExonList().size(); j++) {
                currTranscript.getExonList().get(j).setPos(j);
            }
        }
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

    public HashSet<String> cutCombRegVec(int x1, int x2, AlignmentBlock firstBlock, AlignmentBlock lastBlock) {
        if (combinedExonList == null) {
            initCombinedExonsList();
        }

        HashSet<String> cutRegions = new HashSet<>();
        for (int i = 0; i < combinedExonList.size(); i++) {
            Exon exon = combinedExonList.get(i);
            int exonStart = exon.getStart();
            int exonStop = exon.getStop();

            // #---|-----|------#
            //     x1----x2  → completely contained → add x1-x2 to set
            if (x1 >= exon.getStart() && x1 <= exon.getStop() && x2 >= exon.getStart() && x2 <= exon.getStop()) {
                cutRegions.add(x1 + "-" + x2);
                break;
            }
            // #---------|----|-#
            //           x1---#  → x1 contained → cut
            else if (x1 >= exon.getStart() && x1 <= exon.getStop()) {
                int s = Math.min(exon.getStop(), firstBlock.getReferenceStart() + firstBlock.getLength());
                cutRegions.add(x1 + "-" + s);
            }
            // #----------------#
            // |----------------|
            else if (x1 <= exon.getStart() && x2 >= exon.getStop()) {
                cutRegions.add((exon.getStart()) + "-" + (exon.getStop()));
            }
            // #-------|--------#
            //   X-----x2  → x2 contained → cut
            else if (x2 < exon.getStop()) {
                int s = Math.max(exon.getStart(), lastBlock.getReferenceStart());
                cutRegions.add(s + "-" + x2);
                break;
            }
        }
        return cutRegions;
    }
}

