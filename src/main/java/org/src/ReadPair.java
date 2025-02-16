package org.src;

import augmentedTree.IntervalTree;
import com.sun.source.tree.Tree;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMRecord;

import java.util.*;

public class ReadPair {
    private final SAMRecord fwRecord;
    private final String[] MMKEYWORDS = {"NM", "nM", "XM"};
    private final SAMRecord rwRecord;
    private final Boolean frstrand;
    private final int alignmentStart;
    private final int alignmentEnd;
    private final String chr;
    private final ArrayList<Gene> containingGenes = new ArrayList<>();
    private TreeSet<Region> regionVecFw = new TreeSet<>(Comparator.comparingInt(Region::getStart).thenComparingInt(Region::getStop));
    private TreeSet<Region> regionVecRw = new TreeSet<>(Comparator.comparingInt(Region::getStart).thenComparingInt(Region::getStop));
    private int transcriptomicCount = 0;
    private int mergedCount = 0;
    private int intronCount = 0;
    private int gCount;
    private TreeSet<Region> meltedBlocks;

    public ReadPair(SAMRecord fw, SAMRecord rw, Boolean frstrand) {
        this.fwRecord = fw;
        this.rwRecord = rw;
        this.frstrand = frstrand;
        this.alignmentStart = Math.min(fw.getAlignmentStart(), rw.getAlignmentStart());
        this.alignmentEnd = Math.max(fw.getAlignmentEnd(), rw.getAlignmentEnd());
        this.chr = fw.getReferenceName();
        melt();
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
        // this set will contain all introns
        HashSet<String> iRegions = new HashSet<>();

        // append all implied introns in overlap of fw to iFwRegions
        extractIntrons(overlapStart, overlapEnd, iFwRegions, iRegions, fwRecord);
        // append all implied introns in overlap of rw to iRwRegions
        extractIntrons(overlapStart, overlapEnd, iRwRegions, iRegions, rwRecord);

        // this is a really weird edge case
        if (overlapEnd - overlapStart == -1) {
            return iRegions.size();
        }

        // implied intron missing
        if (iRwRegions.size() != iFwRegions.size()) {
            return -1;
        }

        // if implied introns match → return iRegions
        if (iRwRegions.containsAll(iFwRegions)) {
            return iRegions.size();
        }
        // they dont → inconsistent
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
                .get(this.chr)
                .get(this.frstrand)
                .getIntervalsSpanning(this.alignmentStart, this.alignmentEnd, new ArrayList<>());

        if (!cgenes.isEmpty()) {
            // add genes for later annotation
            containingGenes.addAll(cgenes);
        }
        return cgenes.size();
    }

    public int getigenes(Genome genome) {
        ArrayList<Gene> igenes;
        igenes = genome.getIntervalTreeMap()
                .get(this.chr)
                .get(this.frstrand)
                .getIntervalsSpannedBy(this.alignmentStart, this.alignmentEnd, new ArrayList<>());

        return igenes.size();
    }

    public int getgdist(Genome genome) {
        ArrayList<Gene> leftNeighbors;
        ArrayList<Gene> rightNeighbors;

        leftNeighbors = genome.getIntervalTreeMap()
                .get(this.chr)
                .get(this.frstrand)
                .getIntervalsLeftNeighbor(this.alignmentStart, this.alignmentEnd, new ArrayList<>());
        rightNeighbors = genome.getIntervalTreeMap()
                .get(this.chr)
                .get(this.frstrand)
                .getIntervalsRightNeighbor(this.alignmentStart, this.alignmentEnd, new ArrayList<>());

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
        StringBuilder annotationTranscriptomic = new StringBuilder();
        StringBuilder annotationMerged = new StringBuilder();
        StringBuilder annotationIntronic = new StringBuilder();

        // go through all genes
        for (Gene gene : containingGenes) {
            // check transcriptopmic
            String transcriptomic = findTranscriptomic(gene);
            // handle result
            if (transcriptomic != null) {
                if (!annotationTranscriptomic.isEmpty()) {
                    annotationTranscriptomic.append("|" + transcriptomic);
                    this.transcriptomicCount += 1;
                    continue;
                } else {
                    annotationTranscriptomic.append(transcriptomic);
                    this.transcriptomicCount += 1;
                    continue;
                }
            }


            // check merged
            String merged = findMerged(gene);
            // handle result
            if (merged != null) {
                if (!annotationMerged.isEmpty()) {
                    annotationMerged.append("|" + merged);
                    this.mergedCount += 1;
                    continue;
                } else {
                    annotationMerged.append(merged);
                    this.mergedCount += 1;
                    continue;
                }
            }

            // if not transcriptomic / merged → intronic
            // handle intronic
            if (!annotationIntronic.isEmpty()) {
                annotationIntronic.append("|" + gene.getGeneId() + "," + gene.getBioType() + ":INTRON");
                this.intronCount += 1;
            } else {
                annotationIntronic.append(gene.getGeneId() + "," + gene.getBioType() + ":INTRON");
                this.intronCount += 1;
            }
        }

        // priority → return in order transcriptomic > merged > intronic
        // gCount is updated correspondingly and later written to result
        if (transcriptomicCount != 0) {
            this.gCount = transcriptomicCount;
            return annotationTranscriptomic.toString();
        } else if (mergedCount != 0) {
            this.gCount = mergedCount;
            return annotationMerged.toString();
        } else {
            this.gCount = intronCount;
            return annotationIntronic.toString();
        }
    }

    public String findTranscriptomic(Gene gene) {
        boolean foundTranscript = false;
        StringBuilder transcriptomicSb = new StringBuilder(gene.getGeneId() + "," + gene.getBioType() + ":");
        for (Transcript transcript : gene.getTranscriptList()) {
            // cut fwx1 fwx2 from transcript exons
            TreeSet<Region> cutFwRegions = transcript.cutSet(fwRecord.getAlignmentStart(), fwRecord.getAlignmentEnd());

            // skip if there are no regions in x1-x2
            if (cutFwRegions.isEmpty()) {
                continue;
            }

            // only if the above passes check rw
            if (cutFwRegions.equals(regionVecFw)) {
                // cut rwx1 rwx2 from transcript exons
                TreeSet<Region> cutRwRegions = transcript.cutSet(rwRecord.getAlignmentStart(), rwRecord.getAlignmentEnd());

                // skip if there are no regions in x1-x2
                if (cutRwRegions.isEmpty()) {
                    continue;
                }

                if (cutRwRegions.equals(regionVecRw)) {
                    if (foundTranscript) {
                        transcriptomicSb.append("," + transcript.getTranscriptId());
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
        TreeSet<Region> mergedTranscriptome = gene.getMeltedRegions();
        TreeSet<Region> mergedRead = this.meltedBlocks;

        int count = 0;
        // brute force this one
        for (Region region : mergedRead) {
            for (Region transcriptRegion : mergedTranscriptome) {
                if (transcriptRegion.contains(region)) {
                    count++;
                    break;
                }
            }
        }

        if (!(count == mergedRead.size())) {
            return null;
        }
        return gene.getGeneId() + "," + gene.getBioType() + ":MERGED";
    }

    public int getgCount() {
        return gCount;
    }

    public boolean isAntisense(Genome genome) {
        ArrayList<Gene> cgenes;
        // negate frstrand
        cgenes = genome.getIntervalTreeMap()
                .get(chr)
                .get(!frstrand)
                .getIntervalsSpanning(alignmentStart, alignmentEnd, new ArrayList<>());

        return !(cgenes.isEmpty());
    }

    public TreeSet<Region> getMeltedBlocks() {
        return this.meltedBlocks;
    }

    public void melt() {
        // melt all regions fw, rw, and merged)
        ArrayList<AlignmentBlock> allBlocks = new ArrayList<>();
        ArrayList<AlignmentBlock> fwBlocks = new ArrayList<>();
        ArrayList<AlignmentBlock> rwBlocks = new ArrayList<>();
        fwBlocks.addAll(fwRecord.getAlignmentBlocks());
        rwBlocks.addAll(rwRecord.getAlignmentBlocks());
        allBlocks.addAll(fwBlocks);
        allBlocks.addAll(rwBlocks);
        this.regionVecFw = meltRegion(fwBlocks);
        this.regionVecRw = meltRegion(rwBlocks);
        this.meltedBlocks = meltRegion(allBlocks);
    }

    public TreeSet<Region> meltRegion(ArrayList<AlignmentBlock> blocks) {
        // for a given list of alignmentblocks, melt into TreeSet
        TreeSet<Region> melted = new TreeSet<>(
                Comparator.comparingInt(Region::getStart)
                        .thenComparingInt(Region::getStop)
        );
        if (blocks.size() == 1) {
            melted.add(new Region(blocks.getFirst().getReferenceStart(), blocks.getFirst().getReferenceStart() + blocks.getFirst().getLength() - 1));
            return melted;
        }
        Collections.sort(blocks, Comparator.comparingInt(AlignmentBlock::getReferenceStart));
        AlignmentBlock first = blocks.getFirst();
        Region current = new Region(first.getReferenceStart(), first.getReferenceStart() + first.getLength() - 1);

        for (int i = 1; i < blocks.size(); i++) {
            AlignmentBlock block = blocks.get(i);
            if (block.getReferenceStart() <= current.getStop() + 1) {
                current.setStop(Math.max(current.getStop(), block.getReferenceStart() + block.getLength() - 1));
            } else {
                melted.add(current);
                current = new Region(block.getReferenceStart(), block.getReferenceStart() + block.getLength() - 1);
            }
        }

        // add last block
        melted.add(current);
        return melted;
    }

    public Boolean getFrstrand() {
        return frstrand;
    }
}