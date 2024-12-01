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
            ReadPair pair;

            if (frstrand == null) {
                if (mate.getFirstOfPairFlag()) {
                    pair = new ReadPair(mate, current, null);
                } else {
                    pair = new ReadPair(current, mate, null);
                }
            }
            // fr +
            else if (frstrand) {
                // mate first
                if (mate.getFirstOfPairFlag()) {
                    // mate -
                    if (mate.getReadNegativeStrandFlag()) {
                        pair = new ReadPair(mate, current, false);
                    }
                    // mate +
                    else {
                        pair = new ReadPair(mate, current, true);
                    }
                } else {
                    // curr -
                    if (current.getReadNegativeStrandFlag()) {
                        pair = new ReadPair(current, mate, false);
                    }
                    // curr +
                    else {
                        pair = new ReadPair(current, mate, true);
                    }
                }
            }
            // fr -
            else {
                if (mate.getFirstOfPairFlag()) {
                    // mate -
                    if (mate.getReadNegativeStrandFlag()) {
                        pair = new ReadPair(mate, current, false);
                    }
                    // mate +
                    else {
                        pair = new ReadPair(mate, current, true);
                    }
                } else {
                    // curr -
                    if (current.getReadNegativeStrandFlag()) {
                        pair = new ReadPair(current, mate, false);
                    }
                    // curr +
                    else {
                        pair = new ReadPair(current, mate, true);
                    }
                }
            }


            StringBuilder sb = new StringBuilder(current.getReadName());

            int igenes = pair.getigenes(genome);
            int cgenes = pair.getcgenes(genome);
            int gdist = 0;

            if (cgenes == 0) {
                if (igenes > 0) {
                    continue;
                }
                gdist = pair.getgdist(genome);
                if (gdist == -1) {
                    continue;
                }
            }


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
//                sb.append("\tgcount:0" + "\tgdist:" + gdist);
            } else {
//                sb.append("\tgcount:" + cgenes);
            }
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
}
