package picard.util;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.Usage;

import java.io.File;
import java.io.IOException;

/**
 * @author nhomer
 */
public class BedToIntervalList extends CommandLineProgram {

    @Usage
    public final String USAGE = getStandardUsagePreamble() + "Converts a BED file to an Picard Interval List.";

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input BED file")
    public File INPUT;

    @Option(shortName = StandardOptionDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, doc = "The sequence dictionary")
    public File SEQUENCE_DICTIONARY;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output Picard Interval List")
    public File OUTPUT;

    final Log LOG = Log.getInstance(getClass());

    // Stock main method
    public static void main(final String[] args) {
        new BedToIntervalList().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(SEQUENCE_DICTIONARY);
        IOUtil.assertFileIsWritable(OUTPUT);
        try {
            final SamReader samReader = SamReaderFactory.makeDefault().open(SEQUENCE_DICTIONARY);
            final SAMFileHeader header = samReader.getFileHeader();
            final IntervalList intervalList = new IntervalList(header);
            CloserUtil.close(samReader);

            final FeatureReader<BEDFeature> bedReader = AbstractFeatureReader.getFeatureReader(INPUT.getAbsolutePath(), new BEDCodec(), false);
            final CloseableTribbleIterator<BEDFeature> iterator = bedReader.iterator();
            final ProgressLogger progressLogger = new ProgressLogger(LOG, (int) 1e6);

            while (iterator.hasNext()) {
                final BEDFeature bedFeature = iterator.next();
                final String sequenceName = bedFeature.getChr();
                final int start = bedFeature.getStart(); // 1-based
                final int end = bedFeature.getEnd(); // 1-based
                final SAMSequenceRecord sequenceRecord = header.getSequenceDictionary().getSequence(sequenceName);

                // Do some validation
                if (null == sequenceRecord) {
                    throw new PicardException(String.format("Sequence '%s' was not found in the sequence dictionary", sequenceName));
                }
                else if (start < 1) {
                    throw new PicardException(String.format("Start on sequence '%s' was less than one: %d", sequenceName, start));
                }
                else if (sequenceRecord.getSequenceLength() < start) {
                    throw new PicardException(String.format("Start on sequence '%s' was past the end: %d < %d", sequenceName, sequenceRecord.getSequenceLength(), start));
                }
                else if (end < 1) {
                    throw new PicardException(String.format("End on sequence '%s' was less than one: %d", sequenceName, end));
                }
                else if (sequenceRecord.getSequenceLength() < end) {
                    throw new PicardException(String.format("End on sequence '%s' was past the end: %d < %d", sequenceName, sequenceRecord.getSequenceLength(), end));
                }
                else if (end <= start) {
                    throw new PicardException(String.format("On sequence '%s', end <= start: %d <= %d", sequenceName, end, start));
                }

                final Interval interval = new Interval(sequenceName, start, end);
                intervalList.add(interval);

                progressLogger.record(sequenceName, start);
            }
            CloserUtil.close(bedReader);

            // Sort and write the output
            intervalList.uniqued().write(OUTPUT);

        } catch (final IOException e) {
            throw new RuntimeException(e);
        }

        return 0;
    }
}
