package org.src;

import augmentedTree.IntervalTree;
import net.sf.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;

public class ReadPair {
    private SAMRecord fwRecord;
    private String[] MMKEYWORDS = {"NM", "nM", "XM"};
    private SAMRecord rwRecord;

    private Boolean frstrand;
    private int alignmentStart;
    private int alignmentEnd;
    private String chr;

    public ReadPair(SAMRecord fw, SAMRecord rw, Boolean frstrand) {
        this.fwRecord = fw;
        this.rwRecord = rw;
        this.frstrand = frstrand;
        this.alignmentStart = Math.min(fw.getAlignmentStart(), rw.getAlignmentStart());
        this.alignmentEnd = Math.max(fw.getAlignmentEnd(), rw.getAlignmentEnd());
        this.chr = fw.getReferenceName();
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

    public void extractIntrons(int overlapStart, int overlapEnd, HashSet<String> iRwRegions, HashSet<String> iRegions, SAMRecord rwRecord) {
        // basically extracts introns and adds them to corresponding sets
        for (int i = 0; i < rwRecord.getAlignmentBlocks().size() - 1; i++) {
            int iStart = rwRecord.getAlignmentBlocks().get(i).getReferenceStart() + rwRecord.getAlignmentBlocks().get(i).getLength();
            int iEnd = rwRecord.getAlignmentBlocks().get(i + 1).getReferenceStart();
            if (iStart - iEnd == 0) {
                continue;
            }
            iRegions.add(iStart + "-" + iEnd);
            if ((iStart >= overlapStart && iStart <= overlapEnd) || (iEnd >= overlapStart && iEnd <= overlapEnd)) {
                iRwRegions.add(iStart + "-" + iEnd);
            }
            if (iStart <= overlapStart && iEnd >= overlapEnd) {
                iRwRegions.add(iStart + "-" + iEnd);
            }
        }
    }


    public int getcgenes(Genome genome) {
        ArrayList<Gene> cgenes;

        cgenes = genome.getIntervalTreeMap()
                .get(chr)
                .get(frstrand)
                .getIntervalsSpanning(alignmentStart, alignmentEnd, new ArrayList<>());
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
            if(distance > 0) {
                minDistance = Math.min(minDistance, distance);
            }
        }

        if (minDistance - 1 != 0) {
            return minDistance - 1;
        }

        return -1;
    }
}


//IntervalTree<Intron> intronTree = new IntervalTree<>();
//HashSet<String> iRegions = new HashSet<>();
//HashSet<Integer> iStartsFw = new HashSet<>();
//HashSet<Integer> iEndsFw = new HashSet<>();
//HashSet<Integer> iStartsRw = new HashSet<>();
//HashSet<Integer> iEndsRw = new HashSet<>();
//        for (int i = 0; i < fwRecord.getAlignmentBlocks().size() - 1; i++) {
//int iStart = fwRecord.getAlignmentBlocks().get(i).getReferenceStart() + fwRecord.getAlignmentBlocks().get(i).getLength();
//int iEnd = fwRecord.getAlignmentBlocks().get(i + 1).getReferenceStart();
//Intron intron = new Intron(iStart, iEnd, true);
//            intronTree.add(intron);
//            iRegions.add(iStart + "-" + iEnd);
//            iStartsFw.add(iStart);
//            iEndsFw.add(iEnd);
//        }
//
//                for (int i = 0; i < rwRecord.getAlignmentBlocks().size() - 1; i++) {
//int iStart = rwRecord.getAlignmentBlocks().get(i).getReferenceStart() + rwRecord.getAlignmentBlocks().get(i).getLength();
//int iEnd = rwRecord.getAlignmentBlocks().get(i + 1).getReferenceStart();
//Intron intron = new Intron(iStart, iEnd, false);
//            intronTree.add(intron);
//            iRegions.add(iStart + "-" + iEnd);
//            iStartsRw.add(iStart);
//            iEndsRw.add(iEnd);
//        }
//
//int overlapStart = Math.max(rwRecord.getAlignmentStart(), fwRecord.getAlignmentStart()) + 1;
//int overlapEnd = Math.min(rwRecord.getAlignmentEnd(), fwRecord.getAlignmentEnd());
//
//ArrayList<Intron> overlapIs = intronTree.getIntervalsSpannedBy(overlapStart, overlapEnd, new ArrayList<>());
//
//// I think I need an extra case
//
//// case: no introns in overlap
//        if (overlapIs.isEmpty()) {
//        // there might be introns peeking into overlap
//        for (Intron intron : intronTree) {
//        if ((intron.getStart() >= overlapStart && intron.getStart() <= overlapEnd) || (intron.getStop() <= overlapEnd && intron.getStop() >= overlapStart)) {
//        return -1;
//        }
//        }
//
//        // check for remaining introns if they are peeking into overlap
//        return iRegions.size();
//        }
//                // case: one read has intron in overlap, other one doesnt
//                else if (overlapIs.size() == 1) {
//        return -1;
//        }
//        // case check if all introns have same start / end
//        // meaning: each inton has to have an equivalent partner in the other read
//        else {
//        for (int i = 0; i < overlapIs.size(); i++) {
//Intron currentI = overlapIs.get(i);
//                if (currentI.isFw()) {
//        if (!(iEndsRw.contains(currentI.getStop()) && iStartsRw.contains(currentI.getStart()))) {
//        return -1;
//        }
//        } else {
//        if (!(iEndsFw.contains(currentI.getStop()) && iStartsFw.contains(currentI.getStart()))) {
//        return -1;
//        }
//        }
//        }
//        }
//
//        return iRegions.size();
