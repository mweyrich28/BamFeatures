package org.src;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;


import java.io.IOException;

import static net.sourceforge.argparse4j.impl.Arguments.storeTrue;

public class Main {
    public static void main(String[] args) throws IOException {
        ArgumentParser parser = ArgumentParsers.newFor("ExonSkipping").build().defaultHelp(true).description("Usage:\n\t-gtf <path-to-gtf>\n\t-out <path-to-out-tsv>");
        try {
            parser.addArgument("-gtf").required(true).help("Path to Gene Transfer Format File.");
            parser.addArgument("-bam").required(true).help("Path to Bam File.");
            parser.addArgument("-o").required(true).help("Specify Output File Name.");
            parser.addArgument("-frstrand").help("Specify Strandness of Experiment");

            Namespace ns = parser.parseArgs(args);
            String gtfPath = ns.getString("gtf");
            String out = ns.getString("o");
            String bamPath = ns.getString("bam");
            String st = ns.get("frstrand");

            Boolean strandness = null;
            if (st.equals("true")) {
                strandness = true;
            } else if (st.equals("false")) {
                strandness = false;
            }

//            String pathToBAM = "./BamFeatures/h.sn.1.bam";
//            String pathToGTF = "./BamFeatures/Homo_sapiens.GRCh37.75.gtf";
            BamFeatures bam = new BamFeatures(bamPath, gtfPath, strandness);
            bam.processBAM(strandness, out);
        } catch (ArgumentParserException e) {
            parser.printHelp();
        }
//----------------------------------------------------------------------------------
//        String pathToBAM = "./BamFeatures/y.ns.2.bam";
//        String pathToGTF = "./BamFeatures/Saccharomyces_cerevisiae.R64-1-1.75.gtf";
//        Boolean strandness = null;
//        BamFeatures bam = new BamFeatures(pathToBAM, pathToGTF, strandness);
//----------------------------------------------------------------------------------
//        String pathToBAM = "./BamFeatures/h.sp.3.bam";
//        String pathToGTF = "./BamFeatures/Homo_sapiens.GRCh37.75.gtf";
//        Boolean strandness = true;
//        BamFeatures bam = new BamFeatures(pathToBAM, pathToGTF, strandness);
    }
}