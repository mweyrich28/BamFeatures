package org.src;
import augmentedTree.IntervalTree;
import com.sun.source.tree.Tree;
import net.sf.samtools.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.TreeSet;

public class BamFeatures {

    private final SAMFileReader samReader;
    private final Genome genome;

    public BamFeatures(String pathToBAM, String pathToGTF, Boolean strandness) throws IOException {
        this.genome = new Genome();
        this.genome.readGTF(pathToGTF, strandness);
        this.samReader = new SAMFileReader(new File(pathToBAM), false);
        this.samReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
    }

    public void processBAM(Boolean frstrand, String outPath, boolean lengths) throws IOException {
        HashMap<String, SAMRecord> seenEntries = new HashMap<>();
        HashMap<Boolean, HashMap<TreeSet<Region>, Integer>> pcrIndex = new HashMap<>();
        Iterator<SAMRecord> it = samReader.iterator();
        BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(outPath));
        // track chromosome
        String currentChr = null;
        while (it.hasNext()) {
            SAMRecord current = it.next();

            if (currentChr == null) {
                currentChr = current.getReferenceName();
            } else if (!currentChr.equals(current.getReferenceName())) {
                // clear seen
                seenEntries.clear();
                // clear pcr index
                pcrIndex.clear();
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
            // get mate
            SAMRecord mate = seenEntries.get(current.getReadName());
            // init readpair with its correct strandness
            ReadPair pair = determineReadPair(mate, current, frstrand);

            // append read id to sb
            StringBuilder sb = new StringBuilder(current.getReadName());
            int cgenes = pair.getcgenes(genome);
            int gdist = 0;

            if (cgenes == 0) {
                int igenes = pair.getigenes(genome);
                if (igenes > 0) {
                    continue;
                }
                gdist = pair.getgdist(genome);
            }

            int nsplit = pair.getNsplit();
            if (nsplit == -1) {
                sb.append("\tsplit-inconsistent:true");
                bufferedWriter.write(sb + "\n");
                continue;
            }
            // update count based on annotation
            String annotation = pair.annotateRegion();
            cgenes = pair.getgCount();

            int mm = pair.getmm();
            int clipping = pair.getclipping();

            sb.append("\tmm:" + mm);
            sb.append("\tclipping:" + clipping);
            sb.append("\tnsplit:" + nsplit);

            if (cgenes == 0) {
                sb.append("\tgcount:0" + "\tgdist:" + gdist);
                if (frstrand != null) {
                    sb.append("\tantisense:" + pair.isAntisense(genome));
                } else {
                    sb.append("\tantisense:false");
                }
                // get antsense
            } else {
                sb.append("\tgcount:" + cgenes);
                sb.append("\t" + annotation);
            }
            TreeSet<Region> regions = pair.getMeltedBlocks();
            if (!pcrIndex.containsKey(pair.getFrstrand())) {
                pcrIndex.put(pair.getFrstrand(), new HashMap<>());
            }
            HashMap<TreeSet<Region>, Integer> map = pcrIndex.get(pair.getFrstrand());
            if (map.containsKey(regions)) {
                int count = map.get(regions);
                count++;
                map.put(regions, count);
                sb.append("\tpcrindex: " + count);
            } else {
                map.put(regions, 0);
                sb.append("\tpcrindex: " + 0);
            }

            bufferedWriter.write(sb + "\n");
        }
        bufferedWriter.flush();
        bufferedWriter.close();
        if(lengths) {
            BufferedWriter bf = new BufferedWriter(new FileWriter(outPath +"_gene_lengths.tsv"));
            bf.write("gene\tlength\tchr");
            for(HashMap<Boolean, IntervalTree<Gene>> map : genome.getIntervalTreeMap().values()) {
                for(IntervalTree<Gene> tree: map.values()) {
                    for(Gene gene : tree.descendingSet()) {
                        bf.write("\n" + gene.getGeneId() + "\t" + gene.getMeltedLength() + "\t" + gene.getChr());
                    }
                }
            }
            bf.flush();
            bf.close();
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
                // curr -
                if (current.getReadNegativeStrandFlag()) {
                    return new ReadPair(current, mate, false);
                }
                // curr +
                else {
                    return new ReadPair(current, mate, true);
                }
            } else {
                // mate -
                if (mate.getReadNegativeStrandFlag()) {
                    return new ReadPair(mate, current, false);
                }
                // mate +
                else {
                    return new ReadPair(mate, current, true);
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
