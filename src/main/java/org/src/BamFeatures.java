package org.src;
import net.sf.samtools.*;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class BamFeatures {

    private final SAMFileReader samReader;
    private final Iterator<SAMRecord> it;
    private final Genome genome;

    public BamFeatures(String pathToBAM, String pathToGTF, boolean strandness) throws IOException {
        this.genome = new Genome();
        this.genome.readGTF(pathToGTF);
        this.samReader = new SAMFileReader(new File(pathToBAM), false);
        this.samReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        this.it = samReader.iterator();
        processBAM();
    }

    public void processBAM(){
        HashMap<String, SAMRecord> seenEntries = new HashMap<String, SAMRecord>();
        // track chromosome
        String currentChr = null;
        int count = 0;
        while (it.hasNext()) {
            SAMRecord current = it.next();
            if (currentChr == null) {
                currentChr = current.getReferenceName();
            } else if (!currentChr.equals(current.getReferenceName())) {
                // clear seen
               seenEntries.clear();
            }

            // ignore based on read attributes
            boolean isNotPrimary = current.getNotPrimaryAlignmentFlag();
            boolean isMateMapped = current.getMateUnmappedFlag();
            boolean isNotMapped = current.getReadUnmappedFlag();
            boolean sameStrandChr = current.getReferenceName().equals(current.getMateReferenceName()) && current.getReadNegativeStrandFlag() == (current.getMateNegativeStrandFlag());
            boolean paired = current.getReadPairedFlag();
            if (isNotPrimary|| isNotMapped ||isMateMapped || sameStrandChr || !paired) {
                continue;
            }

            // check if gene is contained inbetween entries
            ArrayList<Gene> genesInInterval = genome.getIntervalTreeMap().get(current.getReferenceName()).get(current.getMateNegativeStrandFlag())
                    .getIntervalsSpanning(current.getAlignmentStart(), current.getMateAlignmentStart(), new ArrayList<>());

            ArrayList<Gene> subInterval = genome.getIntervalTreeMap().get(current.getReferenceName()).get(current.getMateNegativeStrandFlag())
                    .getIntervalsSpannedBy(current.getAlignmentStart(), current.getMateAlignmentStart(), new ArrayList<>());
            boolean containsSingleGene = false;
            for (Gene gene : subInterval) {
                if(gene.getStart() > current.getAlignmentStart() && gene.getEnd() < current.getMateAlignmentStart()) {
                    // TODO: FRAGE muss das nur für ein Gen gelten oder für alle?
                    containsSingleGene = true;
                    break;
                }
            }
            if (containsSingleGene) {
                continue;
            }

            // track entries
            if(!seenEntries.containsKey(current.getReadName())) {
               seenEntries.put(current.getReadName(), current);
               continue;
            }

            // at this point we already have the read pair
            SAMRecord mate = seenEntries.get(current.getReadName());
            ReadPair pair = null;

            if(mate.getFirstOfPairFlag()) {
                pair = new ReadPair(mate, current, this.genome);
            } else {
                pair = new ReadPair(current, mate, this.genome);
            }
            System.out.println(mate.getReadName());
            count++;
        }
        System.out.println();
    }
}
