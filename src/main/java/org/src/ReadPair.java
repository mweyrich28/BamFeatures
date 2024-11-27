package org.src;

import net.sf.samtools.SAMRecord;

public class ReadPair {
    private SAMRecord fwRecord;
    private SAMRecord rwRecord;
    private Genome genome;

    public ReadPair(SAMRecord fw, SAMRecord rw, Genome genome) {
       this.fwRecord = fw;
       this.rwRecord = rw;
       this.genome = genome;
    }

}

