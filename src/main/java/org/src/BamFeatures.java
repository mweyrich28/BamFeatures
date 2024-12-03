package org.src;

import net.sf.samtools.*;
import org.w3c.dom.css.RGBColor;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class BamFeatures {

    private final SAMFileReader samReader;
    private final Genome genome;

    public BamFeatures(String pathToBAM, String pathToGTF, Boolean strandness) throws IOException {
        this.genome = new Genome();
        this.genome.readGTF(pathToGTF, strandness);
        this.samReader = new SAMFileReader(new File(pathToBAM), false);
        this.samReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        processBAM(strandness);
    }

    public void processBAM(Boolean frstrand) {
        HashMap<String, SAMRecord> seenEntries = new HashMap<>();
        Iterator<SAMRecord> it = samReader.iterator();


        // track chromosome
        String currentChr = null;
        while (it.hasNext()) {
            SAMRecord current = it.next();

            if (current.getReadName().equals("2")) {
                System.out.println();
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
            ReadPair pair = determineReadPair(mate, current, frstrand);
            StringBuilder sb = new StringBuilder(current.getReadName());

            int igenes = pair.getigenes(genome);
            int cgenes = pair.getcgenes(genome);
            int gdist = 0;

            if (cgenes == 0) {
                if (igenes > 0) {
                    continue;
                }
                gdist = pair.getgdist(genome);
            }

            pair.annotateRegion(genome);

            int nsplit = pair.getNsplit();
            if (nsplit == -1) {
                sb.append("\tsplit-inconsistent:true");
                System.out.println(sb);
                continue;
            }
            int mm = pair.getmm();
            int clipping = pair.getclipping();

            sb.append("\tmm:" + mm);
            sb.append("\tclipping:" + clipping);
            sb.append("\tnsplit:" + nsplit);

            if (cgenes == 0) {
                sb.append("\tgcount:0" + "\tgdist:" + gdist);
            } else {
                sb.append("\tgcount:" + cgenes);
            }
//            System.out.println(current.getReadName() + "\t");
            System.out.println(sb);
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

    public ReadPair determineReadPair(SAMRecord mate, SAMRecord current, Boolean frstrand) {
        // bases on frstrand, getFirstOfPair and getNegativeStrandFlag determine
        // correct readpair configuration
        if (frstrand == null) {
            if (mate.getFirstOfPairFlag()) {
                return new ReadPair(mate, current, null);
            } else {
                return new ReadPair(current, mate, null);
            }
        }
        // fr +
        else if (frstrand) {
            // mate first
            if (mate.getFirstOfPairFlag()) {
                // mate -
                if (mate.getReadNegativeStrandFlag()) {
                    return new ReadPair(mate, current, false);
                }
                // mate +
                else {
                    return new ReadPair(mate, current, true);
                }
            } else {
                // curr -
                if (current.getReadNegativeStrandFlag()) {
                    return new ReadPair(current, mate, false);
                }
                // curr +
                else {
                    return new ReadPair(current, mate, true);
                }
            }
        }
        // fr -
        else {
            if (mate.getFirstOfPairFlag()) {
                // mate -
                if (mate.getReadNegativeStrandFlag()) {
                    return new ReadPair(mate, current, false);
                }
                // mate +
                else {
                    return new ReadPair(mate, current, true);
                }
            } else {
                // curr -
                if (current.getReadNegativeStrandFlag()) {
                    return new ReadPair(current, mate, false);
                }
                // curr +
                else {
                    return new ReadPair(current, mate, true);
                }
            }
        }
    }
}
