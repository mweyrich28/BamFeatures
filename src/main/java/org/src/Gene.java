package org.src;

import augmentedTree.Interval;

import java.util.ArrayList;
import java.util.HashMap;

public class Gene implements Interval {
    private final int start;
    private final int end;
    private final String geneId;
    private final ArrayList<Transcript> transcriptList;
    private final HashMap<String, Transcript> transcriptMap;
    private final String geneName;
    private final String bioType;
    private final String chr;
    private final char strand;

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

    public void addTranscript(Transcript transcript){
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
}
