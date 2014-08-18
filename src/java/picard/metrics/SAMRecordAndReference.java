package picard.metrics;

import htsjdk.samtools.ReadRecord;
import htsjdk.samtools.reference.ReferenceSequence;

public class SAMRecordAndReference {
    private final ReadRecord samRec;
    private final ReferenceSequence refSeq;

    public SAMRecordAndReference(final ReadRecord samRec, final ReferenceSequence refSeq) {
        this.samRec = samRec;
        this.refSeq = refSeq;
    }

    public ReadRecord getSamRecord() {
        return samRec;
    }

    public ReferenceSequence getReferenceSequence() {
        return refSeq;
    }
}