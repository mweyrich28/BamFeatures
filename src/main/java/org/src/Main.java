package org.src;


import java.io.IOException;

public class Main {
    public static void main(String[] args) throws IOException {
        String pathToBAM = "./BamFeatures/h.sn.1.bam";
        String pathToGTF = "./BamFeatures/Homo_sapiens.GRCh37.75.gtf";
//        String pathToBAM = "./BamFeatures/y.ns.2.bam";
//        String pathToGTF = "./BamFeatures/Saccharomyces_cerevisiae.R64-1-1.75.gtf";
        Boolean strandness = null;
        BamFeatures bam = new BamFeatures(pathToBAM, pathToGTF, strandness);
    }
}