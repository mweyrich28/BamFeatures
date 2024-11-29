package org.src;

import net.sf.samtools.SAMRecord;

import java.util.ArrayList;

public class ReadPair {
    private SAMRecord fwRecord;
    private String[] MMKEYWORDS = {"NM", "nM", "XM"};
    private SAMRecord rwRecord;

    private Boolean frstrand;

    public ReadPair(SAMRecord fw, SAMRecord rw, Boolean frstrand) {
        this.fwRecord = fw;
        this.rwRecord = rw;
        this.frstrand = frstrand;
    }

    public int getmm() {
        int rwmm = 0;
        int fwmm = 0;
//        if (fwRecord.getReadName().equals("10664907")) {
//            System.out.println();
//        }
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
        // no overlap
        if (fwRecord.getAlignmentStart() < rwRecord.getAlignmentStart()) {
            return 0;
        }




        return -1; // split inconsistency
    }
    public String getGenes (Genome genome) {
        ArrayList<Gene> cgenes;
        ArrayList<Gene> igenes;
        ArrayList<Gene> neighboursLeft;
        ArrayList<Gene> neighboursRight;

        // strand unspecific
        if (this.frstrand == null) {
            igenes = genome.getIntervalTreeMap()
                    .get(fwRecord.getReferenceName())
                    .get(fwRecord.getReadNegativeStrandFlag())
                    .getIntervalsSpannedBy(fwRecord.getAlignmentStart(), fwRecord.getMateAlignmentStart(), new ArrayList<>());
            cgenes = genome.getIntervalTreeMap()
                    .get(fwRecord.getReferenceName())
                    .get(fwRecord.getReadNegativeStrandFlag())
                    .getIntervalsSpannedBy(fwRecord.getAlignmentStart(), fwRecord.getMateAlignmentStart(), new ArrayList<>());
            if (!igenes.isEmpty() && cgenes.isEmpty()) {
                return null;
            }
        }
        // + experiment
        else if (frstrand) {

        }
        // - experiment
        else {

        }
        return null;
    }

}

