package org.src;

import net.sf.samtools.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;

public class BamFeatures {

    private final SAMFileReader samReader;
    private final Genome genome;

    public BamFeatures(String pathToBAM, String pathToGTF, Boolean strandness) throws IOException {
        this.genome = new Genome();
        this.genome.readGTF(pathToGTF, strandness);
        this.samReader = new SAMFileReader(new File(pathToBAM), false);
        this.samReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
    }

    public void processBAM(Boolean frstrand, String outPath) throws IOException {
        HashMap<String, SAMRecord> seenEntries = new HashMap<>();
        HashMap<String, Integer> pcrIndex = new HashMap<>();
        Iterator<SAMRecord> it = samReader.iterator();
        BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(new File(outPath)));
        // track chromosome
        String currentChr = null;
        while (it.hasNext()) {
            SAMRecord current = it.next();

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
            if (mate.getReadName().equals("25466612")) {
                System.out.println();
            }

            // merge region vector of reads

            String annotation = pair.annotateRegion();
            // update count based on annotation
            cgenes = pair.getgCount();

            int nsplit = pair.getNsplit();
            if (nsplit == -1) {
                sb.append("\tsplit-inconsistent:true");
                bufferedWriter.write(sb.toString() + "\n");
                continue;
            }
            int mm = pair.getmm();
            int clipping = pair.getclipping();

            sb.append("\tmm:" + mm);
            sb.append("\tclipping:" + clipping);
            sb.append("\tnsplit:" + nsplit);

            if (cgenes == 0) {
                sb.append("\tgcount:0" + "\tgdist:" + gdist);
                if (frstrand != null) {
                    sb.append("\tantisense:" + pair.isAntisense(genome));
                }
                // get antsense
            } else {
                sb.append("\tgcount:" + cgenes);
                sb.append("\t" + annotation);
            }
//            System.out.println(current.getReadName() + "\t");
//            System.out.println(sb);


            String hash = pair.getPcrHash();
            if (pcrIndex.containsKey(hash)) {
                int last = pcrIndex.get(hash);
                pcrIndex.put(hash, last + 1);
                sb.append("\tpcrindex: " + (last + 1));
            } else {
                pcrIndex.put(hash, 0);
                sb.append("\tpcrindex: " + 0);
            }

            bufferedWriter.write(sb.toString() + "\n");
        }
        bufferedWriter.flush();
        bufferedWriter.close();
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
