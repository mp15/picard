package picard.metrics;

import htsjdk.samtools.ReadRecord;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.reference.ReferenceSequence;

/** Defines a MultilevelPerRecordCollector using the argument type of SAMRecord so that this doesn't have to be redefined for each subclass of MultilevelPerRecordCollector */
public abstract class SAMRecordMultiLevelCollector<BEAN extends MetricBase, HKEY extends Comparable> extends MultiLevelCollector<BEAN, HKEY, ReadRecord> {

    @Override
    protected ReadRecord makeArg(ReadRecord samRec, final ReferenceSequence refSeq) {
        return samRec;
    }
}
