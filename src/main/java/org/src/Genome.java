package org.src;

import augmentedTree.IntervalTree;
import org.src.utils.FileUtils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

public class Genome {
    private final HashMap<String, Gene> genes;
    private final HashSet<String> VALID_CHROMOSOMES = new HashSet<>();
    private HashMap<String, HashMap<Boolean, IntervalTree<Gene>>> intervalTreeMap;

    public Genome() {
        this.genes = new HashMap<>();
        this.intervalTreeMap = new HashMap<>();
        initValidChrs();
    }

    public void initValidChrs() {
        this.VALID_CHROMOSOMES.add("1");
        this.VALID_CHROMOSOMES.add("2");
        this.VALID_CHROMOSOMES.add("3");
        this.VALID_CHROMOSOMES.add("4");
        this.VALID_CHROMOSOMES.add("5");
        this.VALID_CHROMOSOMES.add("6");
        this.VALID_CHROMOSOMES.add("7");
        this.VALID_CHROMOSOMES.add("8");
        this.VALID_CHROMOSOMES.add("9");
        this.VALID_CHROMOSOMES.add("10");
        this.VALID_CHROMOSOMES.add("11");
        this.VALID_CHROMOSOMES.add("12");
        this.VALID_CHROMOSOMES.add("13");
        this.VALID_CHROMOSOMES.add("14");
        this.VALID_CHROMOSOMES.add("15");
        this.VALID_CHROMOSOMES.add("16");
        this.VALID_CHROMOSOMES.add("17");
        this.VALID_CHROMOSOMES.add("18");
        this.VALID_CHROMOSOMES.add("19");
        this.VALID_CHROMOSOMES.add("20");
        this.VALID_CHROMOSOMES.add("21");
        this.VALID_CHROMOSOMES.add("22");
        this.VALID_CHROMOSOMES.add("23");
        this.VALID_CHROMOSOMES.add("24");
        this.VALID_CHROMOSOMES.add("25");
        this.VALID_CHROMOSOMES.add("26");
        this.VALID_CHROMOSOMES.add("27");
        this.VALID_CHROMOSOMES.add("28");
        this.VALID_CHROMOSOMES.add("29");
        this.VALID_CHROMOSOMES.add("30");
        this.VALID_CHROMOSOMES.add("31");
        this.VALID_CHROMOSOMES.add("32");
        this.VALID_CHROMOSOMES.add("33");
        this.VALID_CHROMOSOMES.add("34");
        this.VALID_CHROMOSOMES.add("35");
        this.VALID_CHROMOSOMES.add("36");
        this.VALID_CHROMOSOMES.add("37");
        this.VALID_CHROMOSOMES.add("38");
        this.VALID_CHROMOSOMES.add("39");
        this.VALID_CHROMOSOMES.add("40");
        this.VALID_CHROMOSOMES.add("41");
        this.VALID_CHROMOSOMES.add("42");
        this.VALID_CHROMOSOMES.add("43");
        this.VALID_CHROMOSOMES.add("44");
        this.VALID_CHROMOSOMES.add("45");
        this.VALID_CHROMOSOMES.add("46");
        this.VALID_CHROMOSOMES.add("X");
        this.VALID_CHROMOSOMES.add("Y");
        this.VALID_CHROMOSOMES.add("I");
        this.VALID_CHROMOSOMES.add("II");
        this.VALID_CHROMOSOMES.add("III");
        this.VALID_CHROMOSOMES.add("IV");
        this.VALID_CHROMOSOMES.add("IX");
        this.VALID_CHROMOSOMES.add("Mito");
        this.VALID_CHROMOSOMES.add("V");
        this.VALID_CHROMOSOMES.add("VI");
        this.VALID_CHROMOSOMES.add("VII");
        this.VALID_CHROMOSOMES.add("VIII");
        this.VALID_CHROMOSOMES.add("XI");
        this.VALID_CHROMOSOMES.add("XII");
        this.VALID_CHROMOSOMES.add("XIII");
        this.VALID_CHROMOSOMES.add("XIV");
        this.VALID_CHROMOSOMES.add("XV");
        this.VALID_CHROMOSOMES.add("XVI");
        this.VALID_CHROMOSOMES.add("MT");
    }

    public HashMap<String, HashMap<Boolean, IntervalTree<Gene>>> getIntervalTreeMap() {
        return intervalTreeMap;
    }


    public void readGTF(String pathToGtf, Boolean frtrand) throws IOException {
        // sanity check vars
        Gene lastGene = null;
        int exonCounter = 0;

        BufferedReader buff = new BufferedReader(new FileReader(pathToGtf));
        String line;

        while ((line = buff.readLine()) != null) {
            // skip comments
            if (line.charAt(0) == '#') {
                continue;
            }


            // extract main components (line split by \t)
            String[] mainComponents = line.split("\t");
            // skip chromosomes that are not in reference fasta
            if (!VALID_CHROMOSOMES.contains(mainComponents[0])) {
                continue;
            }
            // split attributes again at ";"
            String[] attributeEntries = mainComponents[mainComponents.length - 1].split(";");

            // get newGeneId of current line
            String newGeneId = FileUtils.parseGTFAttributes(attributeEntries, "gene_id");

            // check if we hit a new gene
            if (mainComponents[2].equals("gene")) {
                // update gene and continue with next gtf line
                String bioType = FileUtils.parseGTFAttributes(attributeEntries, "gene_biotype");
                int geneStart = Integer.parseInt(mainComponents[3]);
                int geneEnd = Integer.parseInt(mainComponents[4]);
                String geneName = FileUtils.parseGTFAttributes(attributeEntries, "gene_name");
                String chr = mainComponents[0];
                char strand = mainComponents[6].charAt(0);
                lastGene = new Gene(newGeneId, geneStart, geneEnd, geneName, chr, strand, bioType);
                // false → + strand
                // true  → - strand
                boolean isNegative = strand == '-';

                // Here IntervalTree
                // if frstrand == null → just create one tree
                if (!intervalTreeMap.containsKey(chr)) {
                    intervalTreeMap.put(chr, new HashMap<>());
                }
                if (frtrand == null) {
                    if (!intervalTreeMap.get(chr).containsKey(null)) {
                        intervalTreeMap.get(chr).put(null, new IntervalTree<>());
                    }
                    intervalTreeMap.get(chr).get(null).add(lastGene);

                } else {
                    if (!intervalTreeMap.get(chr).containsKey(isNegative)) {
                        intervalTreeMap.get(chr).put(isNegative, new IntervalTree<>());
                    }
                    intervalTreeMap.get(chr).get(isNegative).add(lastGene);
                }
                continue;
            }


            // did we hit a new transcript
            if (mainComponents[2].equals("transcript")) {

                // only add cds to current transcript
                String transcriptId = FileUtils.parseGTFAttributes(attributeEntries, "transcript_id");

                // add new transcript to current gene
                int transcriptStart = Integer.parseInt(mainComponents[3]);
                int transcriptStop = Integer.parseInt(mainComponents[4]);
                Transcript transcript = new Transcript(transcriptId, mainComponents[6].charAt(0), transcriptStart, transcriptStop);
                lastGene.addTranscript(transcript);

                // reset counter
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
