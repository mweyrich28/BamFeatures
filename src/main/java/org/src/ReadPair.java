package org.src;

import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMRecord;

import java.util.ArrayList;

public class ReadPair {
    private SAMRecord fwRecord;
    private String[] mmKEYWORDS = {"NM", "nM", "XM"};
    private SAMRecord rwRecord;

    private Boolean frstrand;

    public ReadPair(SAMRecord fw, SAMRecord rw, Boolean frstrand) {
        this.fwRecord = fw;
        this.rwRecord = rw;
        this.frstrand = frstrand;
    }

    public int getmm() {
        int mm = 0;
        for (String mmKEYWORD : mmKEYWORDS) {
            if (fwRecord.getAttribute(mmKEYWORD) != null) {
                mm += (Integer) fwRecord.getAttribute(mmKEYWORD);
            }
            if (rwRecord.getAttribute(mmKEYWORD) != null) {
                mm += (Integer) rwRecord.getAttribute(mmKEYWORD);
            }
        }
        return mm;
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

