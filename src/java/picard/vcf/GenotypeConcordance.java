package picard.vcf;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.FormatUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.vcf.VCFFileReader;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.EnumSet;
import java.util.Iterator;
import java.util.List;

/**
 *
 */
public class GenotypeConcordance extends CommandLineProgram {
    @Option(shortName = "TV", doc="The VCF containing the truth sample")
    public File TRUTH_VCF;

    @Option(shortName = "CV", doc="The VCF containing the call sample")
    public File CALL_VCF;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, optional=true, doc="File to output report to, otherwise standard out.")
    public File OUTPUT;

    @Option(shortName = "TS", doc="The name of the truth sample within the truth VCF")
    public String TRUTH_SAMPLE;

    @Option(shortName = "CS", doc="The name of the call sample within the call VCF")
    public String CALL_SAMPLE;

    @Option()
    public List<File> INTERVALS;

    @Option(doc="If true, multiple interval lists will be intersected. If false multiple lists will be unioned.")
    public boolean INTERSECT_INTERVALS = true;

    @Option(doc="Genotypes below this genotype quality will have genotypes classified as LowGq.")
    public int MIN_GQ = 0;

    @Option(doc="Genotypes below this genotype quality will also have genotypes classified as LowDp.")
    public int MIN_DP = 0;

    @Option(doc="If true, output all rows even when count == 0.  When false only output rows with non-zero counts.")
    public boolean OUTPUT_ALL_ROWS = false;

    private final Log log = Log.getInstance(GenotypeConcordance.class);
    private final ProgressLogger progress = new ProgressLogger(log, 10000, "checked", "variants");

    public static void main(final String[] args) {
        new GenotypeConcordance().instanceMainWithExit(args);
    }

    @Override protected int doWork() {
        IOUtil.assertFileIsReadable(TRUTH_VCF);
        IOUtil.assertFileIsReadable(CALL_VCF);
        if (OUTPUT != null) IOUtil.assertFileIsWritable(OUTPUT);
        final boolean usingIntervals = this.INTERVALS != null && this.INTERVALS.size() > 0;
        final VCFFileReader truthReader = new VCFFileReader(TRUTH_VCF, usingIntervals);
        final VCFFileReader callReader = new VCFFileReader(CALL_VCF, usingIntervals);

        // Check that the samples actually exist in the files!
        if (!truthReader.getFileHeader().getSampleNamesInOrder().contains(TRUTH_SAMPLE)) {
            throw new PicardException("File " + TRUTH_VCF.getAbsolutePath() + " does not contain genotypes for sample " + TRUTH_SAMPLE);
        }
        if (!callReader.getFileHeader().getSampleNamesInOrder().contains(CALL_SAMPLE)) {
            throw new PicardException("File " + CALL_VCF.getAbsolutePath() + " does not contain genotypes for sample " + CALL_SAMPLE);
        }

        // TODO: optimization to only use the index if accessing less than < 25% of file?
        // TODO: Add truthSample/callSample detail section.

        // TODO: maybe add optimization if the samples are in the same file??
        // TODO: add option for auto-detect pairs based on same sample name
        // TODO: allow multiple sample-pairs in one pass

        // Build the pair of iterators over the regions of interest
        final Iterator<VariantContext> truthIterator, callIterator;
        if (usingIntervals) {
            log.info("Loading up region lists.");
            IntervalList intervals = null;
            for (final File f : INTERVALS) {
                IOUtil.assertFileIsReadable(f);
                final IntervalList tmp = IntervalList.fromFile(f);

                if (intervals == null)        intervals = tmp;
                else if (INTERSECT_INTERVALS) intervals = IntervalList.intersection(intervals, tmp);
                else intervals =              IntervalList.union(intervals, tmp);
            }

            truthIterator = new MultiIntervalVariantIterator(truthReader, intervals.uniqued());
            callIterator = new MultiIntervalVariantIterator(callReader, intervals.uniqued());
        }
        else {
            truthIterator = truthReader.iterator();
            callIterator = callReader.iterator();
        }

        // Now do the iteration and count things up
        final PairedVariantContextIterator iterator = new PairedVariantContextIterator(truthIterator, callIterator, truthReader.getFileHeader().getSequenceDictionary());
        final ConcordanceResults snpCounter   = new ConcordanceResults(TRUTH_SAMPLE, CALL_SAMPLE);
        final ConcordanceResults indelCounter = new ConcordanceResults(TRUTH_SAMPLE, CALL_SAMPLE);
        log.info("Starting iteration over variants.");

        while (iterator.hasNext()) {
            final VcTuple tuple = iterator.next();
            final State truthState = determineState(tuple.truthVariantContext, TRUTH_SAMPLE, MIN_GQ, MIN_DP);
            final State callState = determineState(tuple.callVariantContext, CALL_SAMPLE, MIN_GQ, MIN_DP);

            if (isSnp(tuple)) snpCounter.add(truthState, callState);
            else indelCounter.add(truthState, callState);

            progress.record(tuple.truthVariantContext == null ? tuple.callVariantContext.getChr() : tuple.truthVariantContext.getChr(), tuple.truthVariantContext == null ? tuple.callVariantContext.getStart() : tuple.truthVariantContext.getStart());
        }

        final PrintStream out;
        if (OUTPUT == null) {
            out = System.out;
        }
        else {
            try { out = new PrintStream(OUTPUT); }
            catch (final IOException ioe) { throw new RuntimeIOException(ioe); }
        }

        final FormatUtil fmt = new FormatUtil();
        out.println("## Summary table:");
        out.println("Event Type\tTruth Sample\tCall Sample\tHet Sens.\tHomVar Sens.\tHet PPV\tHomVar PPV\tVariant Sens.\tVariant PPV");
        out.println("SNP" + "\t" +
                    snpCounter.getTruthSample() + "\t" +
                    snpCounter.getCallSample() + "\t" +
                    fmt.format(snpCounter.hetSensitivity()) + "\t" +
                    fmt.format(snpCounter.homVarSensitivity()) + "\t" +
                    fmt.format(snpCounter.hetPpv()) + "\t" +
                    fmt.format(snpCounter.homVarPpv()) + "\t" +
                    fmt.format(snpCounter.varSensitivity()) + "\t" +
                    fmt.format(snpCounter.varPpv())
        );
        out.println("Indel" + "\t" +
                    indelCounter.getTruthSample() + "\t" +
                    indelCounter.getCallSample() + "\t" +
                    fmt.format(indelCounter.hetSensitivity()) + "\t" +
                    fmt.format(indelCounter.homVarSensitivity()) + "\t" +
                    fmt.format(indelCounter.hetPpv()) + "\t" +
                    fmt.format(indelCounter.homVarPpv()) + "\t" +
                    fmt.format(indelCounter.varSensitivity()) + "\t" +
                    fmt.format(indelCounter.varPpv())
        );
        out.println();

        outputDetailsTable(out, snpCounter, "SNP Detailed Concordance");
        outputDetailsTable(out, indelCounter, "InDel Detailed Concordance");

        return 0;
    }

    /** Outputs the detailed tables for SNP and Indel match categories. */
    private void outputDetailsTable(final PrintStream out, final ConcordanceResults counter, final String description) {
        out.println("## " + description + ":");
        for (final State callState : State.values()) {
            for (final State truthSate : State.values()) {
                    final long count = counter.getCount(callState, truthSate);
                if (count > 0 || OUTPUT_ALL_ROWS) out.println(callState + "\t" + truthSate + "\t" + count);
            }
        }
    }

    /** Determines if the locus is a SNP by querying the first VC first, and only if that is null querying the second VC. */
    final boolean isSnp(final VcTuple tuple) {
        if (tuple.truthVariantContext != null) return tuple.truthVariantContext.isSNP();
        else return tuple.callVariantContext.isSNP();
    }

    /** Determines the classification for a single sample at a single locus. */
    final State determineState(final VariantContext ctx, final String sample, final int minGq, final int minDp) {
        // Site level checks
        if (ctx == null) return State.NoVariant;
        else if (ctx.isFiltered()) return State.FilteredVariant;

        // Genotype level checks
        final Genotype gt = ctx.getGenotype(sample);
        if (gt.isNoCall())           return State.NoCall;
        else if (gt.isFiltered())    return State.FilteredGenotype;
        else if (gt.getGQ() < minGq) return State.LowGq;
        else if (gt.getDP() < minDp) return State.LowDp;
        else if (gt.isHet())         return State.Het;
        else if (gt.isHomRef())      return State.HomRef;
        else if (gt.isHomVar())      return State.HomVar;

        throw new IllegalStateException("Could not classify variant: " + gt);
    }
}

/** Takes a VCFFileReader and an IntervalList and provides a single iterator over all variants in all the intervals. */
class MultiIntervalVariantIterator implements Iterator<VariantContext> {
    private final Iterator<Interval> intervals;
    private final VCFFileReader reader;
    private CloseableIterator<VariantContext> currentIterator;

    MultiIntervalVariantIterator(final VCFFileReader reader, final IntervalList intervals) {
        this.reader = reader;
        this.intervals = intervals.uniqued().iterator();
    }

    /** If the current iterator is null or exhausted, move to the next interval. */
    private void advance() {
        while ((currentIterator == null || !currentIterator.hasNext()) && this.intervals.hasNext()) {
            if (currentIterator != null) currentIterator.close();
            final Interval i = this.intervals.next();
            this.currentIterator = this.reader.query(i.getSequence(), i.getStart(), i.getEnd());
        }
    }

    @Override public boolean hasNext() {
        advance();
        return this.currentIterator.hasNext();
    }

    @Override public VariantContext next() {
        advance();
        return this.currentIterator.next();
    }

    @Override public void remove() {
        throw new UnsupportedOperationException();
    }
}

/**
 * Enum class that provides a total classification for a call, or lack thereof, for a sample at a locus.
 */
enum State {HomRef, Het, HomVar, FilteredVariant, FilteredGenotype, LowGq, LowDp, NoCall, NoVariant}

/** Little class to hold a pair of VariantContexts that are in sync with one another. */
class VcTuple {
    public final VariantContext truthVariantContext;
    public final VariantContext callVariantContext;

    VcTuple(final VariantContext truthVariantContext, final VariantContext callVariantContext) {
        this.truthVariantContext = truthVariantContext;
        this.callVariantContext = callVariantContext;
    }
}

/** Iterator that takes a pair of iterators over VariantContexts and iterates over them in tandem. */
class PairedVariantContextIterator implements Iterator<VcTuple> {
    private final PeekableIterator<VariantContext> truthIterator;
    private final PeekableIterator<VariantContext> callIterator;
    private final VariantContextComparator comparator;

    PairedVariantContextIterator(final Iterator<VariantContext> truthIterator, final Iterator<VariantContext> callIterator, final SAMSequenceDictionary dict) {
        this.truthIterator = new PeekableIterator<VariantContext>(truthIterator);
        this.callIterator = new PeekableIterator<VariantContext>(callIterator);
        this.comparator = new VariantContextComparator(dict);
    }

    @Override
    public boolean hasNext() {
        return this.truthIterator.hasNext() || this.callIterator.hasNext();
    }

    @Override
    public VcTuple next() {
        if (!hasNext()) throw new IllegalStateException("next() called while hasNext() is false.");

        final VariantContext truthIterator = this.truthIterator.hasNext() ? this.truthIterator.peek() : null;
        final VariantContext callIterator = this.callIterator.hasNext() ? this.callIterator.peek() : null;

        // If one or the other is null because there is no next, just return a one-sided tuple
        if (truthIterator == null)       return new VcTuple(null,                  this.callIterator.next());
        else if (callIterator == null)  return new VcTuple(this.truthIterator.next(), null);

        // Otherwise check the ordering and do the right thing
        final int ordering = this.comparator.compare(truthIterator, callIterator);
        if (ordering == 0)     return new VcTuple(this.truthIterator.next(), this.callIterator.next());
        else if (ordering < 0) return new VcTuple(this.truthIterator.next(), null);
        else                   return new VcTuple(null,                  this.callIterator.next());
    }

    @Override public void remove() {
        throw new UnsupportedOperationException();
    }
}

class ConcordanceResults {
    private static class TwoState implements Comparable<TwoState> {
        final State truthState, callState;

        private TwoState(final State truthState, final State callState) {
            this.truthState = truthState;
            this.callState = callState;
        }

        @Override public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            return compareTo((TwoState) o) == 0;
        }

        @Override public int hashCode() {
            int result = truthState.hashCode();
            result = 31 * result + callState.hashCode();
            return result;
        }

        @Override public int compareTo(final TwoState that) {
            int result = this.truthState.compareTo(that.truthState);
            if (result == 0) result = this.callState.compareTo(that.callState);
            return result;
        }
    }

    private final Histogram<TwoState> counts = new Histogram<TwoState>();
    private final String truthSample, callSample;

    ConcordanceResults(final String truthSample, final String callSample) {
        this.truthSample = truthSample;
        this.callSample = callSample;
    }

    void add(final State truthState, final State callState) {
        this.counts.increment(new TwoState(truthState, callState));
    }

    long getCount(final State truthState, final State callState) {
        final Histogram<TwoState>.Bin bin = this.counts.get(new TwoState(truthState, callState));
        if (bin == null) return 0;
        else return (long) bin.getValue();
    }

    public String getCallSample() { return callSample; }
    public String getTruthSample() { return truthSample; }

    double hetSensitivity() {
        return fraction(EnumSet.of(State.Het), EnumSet.of(State.Het),
                        EnumSet.of(State.Het), EnumSet.allOf(State.class));
    }

    double homVarSensitivity() {
        return fraction(EnumSet.of(State.HomVar), EnumSet.of(State.HomVar),
                        EnumSet.of(State.HomVar), EnumSet.allOf(State.class));
    }

    double hetPpv() {
        return fraction(EnumSet.of(State.Het),      EnumSet.of(State.Het),
                        EnumSet.allOf(State.class), EnumSet.of(State.Het));
    }

    double homVarPpv() {
        return fraction(EnumSet.of(State.HomVar),   EnumSet.of(State.HomVar),
                        EnumSet.allOf(State.class), EnumSet.of(State.HomVar));
    }

    double varSensitivity() {
        return fraction(EnumSet.of(State.Het, State.HomVar),   EnumSet.of(State.Het, State.HomVar),
                        EnumSet.of(State.Het, State.HomVar),   EnumSet.allOf(State.class));
    }

    double varPpv() {
        return fraction(EnumSet.of(State.Het, State.HomVar),   EnumSet.of(State.Het, State.HomVar),
                        EnumSet.allOf(State.class),            EnumSet.of(State.Het, State.HomVar));
    }

    private double fraction(final EnumSet<State> lhsNum, final EnumSet<State> rhsNum, final EnumSet<State> lhsDenom, final EnumSet<State> rhsDenom) {
        return sum(lhsNum, rhsNum) / sum(lhsDenom, rhsDenom);
    }

    /** Sums the counts where the first state is contains in truthStateSet and the second state is contained in callStateSet. */
    private double sum(final EnumSet<State> truthStateSet, final EnumSet<State> callStateSet) {
        double result = 0;
        for (final State truthState : truthStateSet) {
            for (final State callState : callStateSet) {
                result += getCount(truthState, callState);
            }
        }

        return result;
    }

}

