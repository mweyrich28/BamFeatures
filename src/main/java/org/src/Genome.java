package org.src;

import augmentedTree.IntervalTree;
import org.src.utils.FileUtils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

public class Genome {
    private final HashMap<String, Gene> genes;
    private HashMap<String, HashMap<Boolean , IntervalTree<Gene>>> intervalTreeMap;
    public Genome() {
        this.genes = new HashMap<>();
        this.intervalTreeMap = new HashMap<>();
    }

    public HashMap<String, HashMap<Boolean, IntervalTree<Gene>>> getIntervalTreeMap() {
        return intervalTreeMap;
    }

    public HashMap<String, Gene> getGenes() {
        return genes;
    }

    public void readGTF(String pathToGtf) throws IOException {
        // sanity check vars
        Gene lastGene = null;
        int exonCounter = 0;

        BufferedReader buff = new BufferedReader(new FileReader(pathToGtf));
        String line;

        while((line = buff.readLine()) != null) {
            // skip comments
            if (line.charAt(0) == '#') {
                continue;
            }


            // extract main components (line split by \t)
            String[] mainComponents = line.split("\t");
            // split attributes again at ";"
            String[] attributeEntries = mainComponents[mainComponents.length - 1].split(";");

            // get newGeneId of current line
            String newGeneId = FileUtils.parseGTFAttributes(attributeEntries, "gene_id");

            // check if we hit a new gene
            if (mainComponents[2].equals("gene")) {
                // update gene and continue with next gtf line
                String bioType = mainComponents[1];
                int geneStart = Integer.parseInt(mainComponents[3]);
                int geneEnd = Integer.parseInt(mainComponents[4]);
                String geneName = FileUtils.parseGTFAttributes(attributeEntries, "gene_name");
                String chr = mainComponents[0];
                char strand = mainComponents[6].charAt(0);
                lastGene = new Gene(newGeneId, geneStart, geneEnd, geneName, chr, strand, bioType);
                // false → + strand
                // true  → - strand
                boolean isNegative = strand == '-';

                // TODO: Keep for now
                genes.put(lastGene.getGeneId(), lastGene);

                // Here IntervalTree ds
                if (!intervalTreeMap.containsKey(chr)) {
                    intervalTreeMap.put(chr, new HashMap<>());
                }
                if (!intervalTreeMap.get(chr).containsKey(isNegative)) {
                    intervalTreeMap.get(chr).put(isNegative, new IntervalTree<>());
                }
                intervalTreeMap.get(chr).get(isNegative).add(lastGene);
                continue;
            }


            // did we hit a new transcript
            if (mainComponents[2].equals("transcript")) {

                // only add cds to current transcript
                String transcriptId = FileUtils.parseGTFAttributes(attributeEntries, "transcript_id");

                // add gene to genome
                genes.put(lastGene.getGeneId(), lastGene);

                // add new transcript to current gene
                int transcriptStart = Integer.parseInt(mainComponents[3]);
                int transcriptStop = Integer.parseInt(mainComponents[4]);
                Transcript transcript = new Transcript(transcriptId, mainComponents[2], transcriptStart, transcriptStop);
                lastGene.addTranscript(transcript);

                // reset
                exonCounter = 0;
            }
            // add exon to last transcript
            else if (mainComponents[2].equals("exon")) {
                int start = Integer.parseInt(mainComponents[3]);
                int end = Integer.parseInt(mainComponents[4]);
                lastGene.getLastTranscript().addExon(
                        start,
                        end,
                        exonCounter
                );
                exonCounter++;
            }
        }
    }
}
