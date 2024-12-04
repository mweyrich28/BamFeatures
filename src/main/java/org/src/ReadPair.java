package org.src;

import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.HashSet;

public class ReadPair {
    private final SAMRecord fwRecord;
    private final String[] MMKEYWORDS = {"NM", "nM", "XM"};
    private final SAMRecord rwRecord;

    private final Boolean frstrand;
    private final int alignmentStart;
    private final int alignmentEnd;
    private final String chr;
    private final ArrayList<Gene> containingGenes = new ArrayList<>();
    private final HashSet<String> regionVecFw = new HashSet<>();
    private final HashSet<String> regionVecRw = new HashSet<>();

    public ReadPair(SAMRecord fw, SAMRecord rw, Boolean frstrand) {
        this.fwRecord = fw;
        this.rwRecord = rw;
        this.frstrand = frstrand;
        this.alignmentStart = Math.min(fw.getAlignmentStart(), rw.getAlignmentStart());
        this.alignmentEnd = Math.max(fw.getAlignmentEnd(), rw.getAlignmentEnd());
        this.chr = fw.getReferenceName();

        // init region vec
        for (AlignmentBlock block : fwRecord.getAlignmentBlocks()) {
            regionVecFw.add(block.getReferenceStart() + "-" + (block.getReferenceStart() + block.getLength() - 1));
        }
        for (AlignmentBlock block : rwRecord.getAlignmentBlocks()) {
            regionVecRw.add(block.getReferenceStart() + "-" + (block.getReferenceStart() + block.getLength() - 1));
        }
    }

    public int getmm() {
        int rwmm = 0;
        int fwmm = 0;
        for (String key : MMKEYWORDS) {
            if (fwRecord.getAttribute(key) != null && (Integer) fwRecord.getAttribute(key) > fwmm) {
                fwmm = (Integer) fwRecord.getAttribute(key);
            }
            if (rwRecord.getAttribute(key) != null && (Integer) rwRecord.getAttribute(key) > rwmm) {
                rwmm = (Integer) rwRecord.getAttribute(key);
            }
        }
        return fwmm + rwmm;
    }

    public int getclipping() {
        return fwRecord.getAlignmentStart() - fwRecord.getUnclippedStart()
                + fwRecord.getUnclippedEnd() - fwRecord.getAlignmentEnd()
                + rwRecord.getAlignmentStart() - rwRecord.getUnclippedStart()
                + rwRecord.getUnclippedEnd() - rwRecord.getAlignmentEnd();
    }

    public int getNsplit() {
        // no splits in records
        if (fwRecord.getAlignmentBlocks().size() == 1 && rwRecord.getAlignmentBlocks().size() == 1) {
            return 0;
        }

        // get overlap
        int overlapStart = Math.max(rwRecord.getAlignmentStart(), fwRecord.getAlignmentStart()) + 1;
        int overlapEnd = Math.min(rwRecord.getAlignmentEnd(), fwRecord.getAlignmentEnd());
        // sometimes swap is needed
        if (overlapStart > overlapEnd) {
            int tmp = overlapEnd;
            overlapEnd = overlapStart - 1;
            overlapStart = tmp + 1;
        }

        HashSet<String> iFwRegions = new HashSet<>();
        HashSet<String> iRwRegions = new HashSet<>();
        HashSet<String> iRegions = new HashSet<>();
        extractIntrons(overlapStart, overlapEnd, iFwRegions, iRegions, fwRecord);

        extractIntrons(overlapStart, overlapEnd, iRwRegions, iRegions, rwRecord);

        // implied intron missing
        if (iRwRegions.size() != iFwRegions.size()) {
            return -1;
        }

        // overlapping introns not matching
        if (iRwRegions.containsAll(iFwRegions)) {
            return iRegions.size();
        }
        return -1;
    }

    public void extractIntrons(int overlapStart, int overlapEnd, HashSet<String> recordRegions, HashSet<String> iRegions, SAMRecord record) {
        // basically extracts introns and adds them to corresponding sets
        for (int i = 0; i < record.getAlignmentBlocks().size() - 1; i++) {
            int iStart = record.getAlignmentBlocks().get(i).getReferenceStart() + record.getAlignmentBlocks().get(i).getLength();
            int iEnd = record.getAlignmentBlocks().get(i + 1).getReferenceStart();
            if (iStart - iEnd == 0) {
                continue;
            }
            iRegions.add(iStart + "-" + iEnd);
            if ((iStart >= overlapStart && iStart <= overlapEnd) || (iEnd >= overlapStart && iEnd <= overlapEnd)) {
                recordRegions.add(iStart + "-" + iEnd);
            }
            if (iStart <= overlapStart && iEnd >= overlapEnd) {
                recordRegions.add(iStart + "-" + iEnd);
            }
        }
    }

    public int getcgenes(Genome genome) {
        ArrayList<Gene> cgenes;
        cgenes = genome.getIntervalTreeMap()
                .get(chr)
                .get(frstrand)
                .getIntervalsSpanning(alignmentStart, alignmentEnd, new ArrayList<>());

        if (!cgenes.isEmpty()) {
            // add genes for later
            containingGenes.addAll(cgenes);
        }
        return cgenes.size();
    }

    public int getigenes(Genome genome) {
        ArrayList<Gene> igenes;
        igenes = genome.getIntervalTreeMap()
                .get(chr)
                .get(frstrand)
                .getIntervalsSpannedBy(alignmentStart, alignmentEnd, new ArrayList<>());
        return igenes.size();
    }

    public int getgdist(Genome genome) {
        ArrayList<Gene> leftNeighbors;
        ArrayList<Gene> rightNeighbors;

        leftNeighbors = genome.getIntervalTreeMap()
                .get(chr)
                .get(frstrand)
                .getIntervalsLeftNeighbor(alignmentStart, alignmentEnd, new ArrayList<>());
        rightNeighbors = genome.getIntervalTreeMap()
                .get(chr)
                .get(frstrand)
                .getIntervalsRightNeighbor(alignmentStart, alignmentEnd, new ArrayList<>());

        // min(read start - gene end  or gene start - read end)
        int minDistance = Integer.MAX_VALUE;
        for (Gene gene : leftNeighbors) {
            int distance = alignmentStart - gene.getEnd(); // read start - gene end
            if (distance > 0) {
                minDistance = Math.min(minDistance, distance);
            }
        }

        for (Gene gene : rightNeighbors) {
            int distance = gene.getStart() - alignmentEnd; // gene start - read end
            if (distance > 0) {
                minDistance = Math.min(minDistance, distance);
            }
        }

        if (minDistance == Integer.MAX_VALUE) {
            return 0;
        }

        return minDistance - 1;
    }

    public String annotateRegion() {
        for (Gene gene : containingGenes) {
            String transcriptomic = findTranscriptomic(gene);
            if (transcriptomic != null) {
                return transcriptomic;
            }


            if (fwRecord.getReadName().equals("361947")) {
                System.out.println();
            }
            String merged = findMerged(gene);
            if (merged != null) {
               return merged;
            }
            return gene.getGeneId() + "," + gene.getBioType() + ":INTRON";
        }
        return null;
    }

    public String findTranscriptomic(Gene gene) {
        boolean foundTranscript = false;
        StringBuilder transcriptomicSb = new StringBuilder(gene.getGeneId() + "," + gene.getBioType() + ":");
        for (Transcript transcript : gene.getTranscriptList()) {
            HashSet<String> trRegionsFw = transcript.cut(fwRecord.getAlignmentStart(), fwRecord.getAlignmentEnd());
            if (regionVecFw.equals(trRegionsFw)) {
                HashSet<String> trRegionsRw = transcript.cut(rwRecord.getAlignmentStart(), rwRecord.getAlignmentEnd());
                if (regionVecRw.equals(trRegionsRw)) {
                    if (foundTranscript) {
                        transcriptomicSb.append("|" + transcript.getTranscriptId());
                    } else {
                        transcriptomicSb.append(transcript.getTranscriptId());
                        foundTranscript = true;
                    }
                }
            }
        }
        if (foundTranscript) {
            return transcriptomicSb.toString();
        } else {
            return null;
        }
    }

    public String findMerged(Gene gene) {
        // check fwRead
        boolean foundMergedFw = false;
        boolean foundMergedRw = false;
        HashSet<String> combRegVec = gene.cutCombRegVec(fwRecord.getAlignmentStart(), fwRecord.getAlignmentEnd(), fwRecord.getAlignmentBlocks().getFirst(), fwRecord.getAlignmentBlocks().getLast());
        // handle edge case where combRegVec.size() == 1 →
        // cR  x1------------------------------------x2
        // fR  x1----------#  #-----#     #----------x2
        if (combRegVec.size() == 1) {
            foundMergedFw = true;
        }
        else if (combRegVec.containsAll(regionVecFw)) {
            foundMergedFw = true;
        }

        if(foundMergedFw) {
            // if there is more than 1 transcript, get combRegVec of Gene
            combRegVec = gene.cutCombRegVec(rwRecord.getAlignmentStart(), rwRecord.getAlignmentEnd(), rwRecord.getAlignmentBlocks().getFirst(), rwRecord.getAlignmentBlocks().getLast());
            if (combRegVec.size() == 1) {
                foundMergedRw = true;
            }
            else if (combRegVec.containsAll(regionVecRw)) {
                foundMergedRw = true;
            }
        }

        if (foundMergedFw && foundMergedRw) {
            return gene.getGeneId() + "," + gene.getBioType() + ":MERGED";
        } else {
            return null;
        }
    }
}

