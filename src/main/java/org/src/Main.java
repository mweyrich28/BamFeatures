package org.src;


import java.io.IOException;

public class Main {
    public static void main(String[] args) throws IOException {
        String pathToBAM = "./BamFeatures/h.sn.1.bam";
        String pathToGTF = "./BamFeatures/Homo_sapiens.GRCh37.75.gtf";
        BamFeatures bam = new BamFeatures(pathToBAM, pathToGTF, false);
    }
}