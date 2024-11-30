package org.src;

import net.sf.samtools.*;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class BamFeatures {

    private final SAMFileReader samReader;
    private final Iterator<SAMRecord> it;
    private final Genome genome;

    public BamFeatures(String pathToBAM, String pathToGTF, Boolean strandness) throws IOException {
        this.genome = new Genome();
        this.genome.readGTF(pathToGTF);
        this.samReader = new SAMFileReader(new File(pathToBAM), false);
        this.samReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        this.it = samReader.iterator();
        processBAM(strandness);
    }

    public void processBAM(Boolean frstrand) {
        HashMap<String, SAMRecord> seenEntries = new HashMap<>();
        // track chromosome
        String currentChr = null;
        while (it.hasNext()) {
            SAMRecord current = it.next();
            if(current.getReadName().equals("59604")) {
//                System.out.println();
            }
            if (currentChr == null) {
                currentChr = current.getReferenceName();
            } else if (!currentChr.equals(current.getReferenceName())) {
                // clear seen
                seenEntries.clear();
                // clear chr in intervaltree as well
                genome.getIntervalTreeMap().get(currentChr).clear();
                // update currChr
                currentChr = current.getReferenceName();
            }

            if (!flagCheck(current)) {
                continue;
            }

            // track entries
            if (!seenEntries.containsKey(current.getReadName())) {
                seenEntries.put(current.getReadName(), current);
                continue;
            }

            // at this point we already have the read pair
            SAMRecord mate = seenEntries.get(current.getReadName());
            ReadPair pair;

            if (mate.getAlignmentStart() < current.getAlignmentStart()) {
                pair = new ReadPair(mate, current, frstrand);
            } else {
                pair = new ReadPair(current, mate, frstrand);
            }

            if (!pair.preGeneCheck(genome)){
               continue;
            }

            int nsplit = pair.getNsplit();
            if (nsplit == -1) {
                System.out.println(current.getReadName() + "\tsplit-inconsistent:true");
                continue;
            }
            int mm = pair.getmm();
            int clipping = pair.getclipping();
            System.out.println(current.getReadName() + "\tmm:" + mm + "\tclipping:" + clipping + "\tnsplit:" + nsplit);
//            System.out.println(current.getReadName() + "\tmm:" + mm + "\tclipping:" + clipping);
//            System.out.println(current.getReadName() + "\tclipping:" + clipping);
//            System.out.println(current.getReadName());
        }
    }
    public boolean flagCheck(SAMRecord record) {
        // ignore based on flags
        boolean isPrimary = !record.getNotPrimaryAlignmentFlag();
        boolean isMateMapped = !record.getMateUnmappedFlag();
        boolean isMapped = !record.getReadUnmappedFlag();
        boolean sameChr = record.getReferenceName().equals(record.getMateReferenceName());
        boolean oppStrand = record.getReadNegativeStrandFlag() != record.getMateNegativeStrandFlag();
        boolean paired = record.getReadPairedFlag();
        return isPrimary && isMapped && isMateMapped && sameChr && oppStrand && paired;
    }

    public boolean preGeneCheck(SAMRecord record) {
        ArrayList<Gene> igenes = genome.getIntervalTreeMap()
                .get(record.getReferenceName())
                .get(record.getReadNegativeStrandFlag())
                .getIntervalsSpannedBy(record.getAlignmentStart(), record.getMateAlignmentStart(), new ArrayList<>());
        ArrayList<Gene> cgenes = genome.getIntervalTreeMap()
                .get(record.getReferenceName())
                .get(record.getReadNegativeStrandFlag())
                .getIntervalsSpanning(record.getAlignmentStart(), record.getMateAlignmentStart(), new ArrayList<>());
        return igenes.isEmpty() || !cgenes.isEmpty();
    }

}
