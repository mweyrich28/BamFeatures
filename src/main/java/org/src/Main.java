package org.src;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;
import java.io.IOException;

import static net.sourceforge.argparse4j.impl.Arguments.storeTrue;

public class Main {
    public static void main(String[] args) throws IOException {
        ArgumentParser parser = ArgumentParsers.newFor("BamFeatures").build().defaultHelp(true).description("Usage:\n\t-gtf <path-to-gtf>\n\t-o <path-to-out.annot>\n\t-bam <path-to-bam>\n\t[-frstrand <true/false>]\n\t-lengths");
        try {
            parser.addArgument("-gtf").required(true).help("Path to Gene Transfer Format File.");
            parser.addArgument("-bam").required(true).help("Path to Bam File.");
            parser.addArgument("-o").required(true).help("Specify Output File Name.");
            parser.addArgument("-frstrand").help("Specify Strandness of Experiment. [-frstrand: true/false]");
            parser.addArgument("-lengths").action(storeTrue()).help("Used to get merged transcript length of genes in ouptut");

            Namespace ns = parser.parseArgs(args);
            String gtfPath = ns.getString("gtf");
            String out = ns.getString("o");
            String bamPath = ns.getString("bam");
            String st = ns.get("frstrand");
            boolean lengths = ns.getBoolean("lengths");

            Boolean strandness = null;
            if (st != null) {
                if (st.equals("true")) {
                    strandness = true;
                } else if (st.equals("false")) {
                    strandness = false;
                }
            }
            BamFeatures bam = new BamFeatures(bamPath, gtfPath, strandness);
            bam.processBAM(strandness, out, lengths);
        } catch (ArgumentParserException e) {
            parser.printHelp();
        }
    }
}