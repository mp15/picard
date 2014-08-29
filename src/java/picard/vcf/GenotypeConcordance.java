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
    @Option(shortName = "V1")
    public File VCF1;

    @Option(shortName = "V2")
    public File VCF2;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, optional=true, doc="File to output report to, otherwise standard out.")
    public File OUTPUT;

    @Option(shortName = "S1")
    public String SAMPLE1;

    @Option(shortName = "S2")
    public String SAMPLE2;

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
        IOUtil.assertFileIsReadable(VCF1);
        IOUtil.assertFileIsReadable(VCF2);
        if (OUTPUT != null) IOUtil.assertFileIsWritable(OUTPUT);
        final boolean usingIntervals = this.INTERVALS != null && this.INTERVALS.size() > 0;
        final VCFFileReader r1 = new VCFFileReader(VCF1, usingIntervals);
        final VCFFileReader r2 = new VCFFileReader(VCF2, usingIntervals);

        // Check that the samples actually exist in the files!
        if (!r1.getFileHeader().getSampleNamesInOrder().contains(SAMPLE1)) {
            throw new PicardException("File " + VCF1.getAbsolutePath() + " does not contain genotypes for sample " + SAMPLE1);
        }
        if (!r2.getFileHeader().getSampleNamesInOrder().contains(SAMPLE2)) {
            throw new PicardException("File " + VCF2.getAbsolutePath() + " does not contain genotypes for sample " + SAMPLE2);
        }

        // TODO: optimization to only use the index if accessing less than < 25% of file?
        // TODO: Add sample1/sample2 detail section.

        // TODO: maybe add optimization if the samples are in the same file??
        // TODO: add option for auto-detect pairs based on same sample name
        // TODO: allow multiple sample-pairs in one pass

        // Build the pair of iterators over the regions of interest
        final Iterator<VariantContext> iterator1, iterator2;
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

            iterator1 = new MultiIntervalVariantIterator(r1, intervals.uniqued());
            iterator2 = new MultiIntervalVariantIterator(r2, intervals.uniqued());
        }
        else {
            iterator1 = r1.iterator();
            iterator2 = r2.iterator();
        }

        // Now do the iteration and count things up
        final PairedVariantContextIterator iterator = new PairedVariantContextIterator(iterator1, iterator2, r1.getFileHeader().getSequenceDictionary());
        final ConcordanceResults snpCounter   = new ConcordanceResults(SAMPLE1, SAMPLE2);
        final ConcordanceResults indelCounter = new ConcordanceResults(SAMPLE1, SAMPLE2);
        log.info("Starting iteration over variants.");

        while (iterator.hasNext()) {
            final VcTuple tuple = iterator.next();
            final State state1 = determineState(tuple.vc1, SAMPLE1, MIN_GQ, MIN_DP);
            final State state2 = determineState(tuple.vc2, SAMPLE2, MIN_GQ, MIN_DP);

            if (isSnp(tuple)) snpCounter.add(state1, state2);
            else indelCounter.add(state1, state2);

            progress.record(tuple.vc1 == null ? tuple.vc2.getChr() : tuple.vc1.getChr(), tuple.vc1 == null ? tuple.vc2.getStart() : tuple.vc1.getStart());
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
        out.println("Event Type\tSample1\tSample2\tHet Sens.\tHomVar Sens.\tHet PPV\tHomVar PPV\tVariant Sens.\tVariant PPV");
        out.println("SNP" + "\t" +
                    snpCounter.getSample1() + "\t" +
                    snpCounter.getSample2() + "\t" +
                    fmt.format(snpCounter.hetSensitivity()) + "\t" +
                    fmt.format(snpCounter.homVarSensitivity()) + "\t" +
                    fmt.format(snpCounter.hetPpv()) + "\t" +
                    fmt.format(snpCounter.homVarPpv()) + "\t" +
                    fmt.format(snpCounter.varSensitivity()) + "\t" +
                    fmt.format(snpCounter.varPpv())
        );
        out.println("Indel" + "\t" +
                    indelCounter.getSample1() + "\t" +
                    indelCounter.getSample2() + "\t" +
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
        for (final State state1 : State.values()) {
            for (final State state2 : State.values()) {
                final long count = counter.getCount(state1, state2);
                if (count > 0 || OUTPUT_ALL_ROWS) out.println(state1 + "\t" + state2 + "\t" + count);
            }
        }
    }

    /** Determines if the locus is a SNP by querying the first VC first, and only if that is null querying the second VC. */
    final boolean isSnp(final VcTuple tuple) {
        if (tuple.vc1 != null) return tuple.vc1.isSNP();
        else return tuple.vc2.isSNP();
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
enum State {HomRef, Het, HomVar, FilteredVariant, FilteredGenotype, LowGq, LowDp, NoCall, NoVariant};

/** Little class to hold a pair of VariantContexts that are in sync with one another. */
class VcTuple {
    public final VariantContext vc1;
    public final VariantContext vc2;

    VcTuple(final VariantContext vc1, final VariantContext vc2) {
        this.vc1 = vc1;
        this.vc2 = vc2;
    }
}

/** Iterator that takes a pair of iterators over VariantContexts and iterates over them in tandem. */
class PairedVariantContextIterator implements Iterator<VcTuple> {
    private final PeekableIterator<VariantContext> iterator1;
    private final PeekableIterator<VariantContext> iterator2;
    private final VariantContextComparator comparator;

    PairedVariantContextIterator(final Iterator<VariantContext> iterator1, final Iterator<VariantContext> iterator2, final SAMSequenceDictionary dict) {
        this.iterator1 = new PeekableIterator<VariantContext>(iterator1);
        this.iterator2 = new PeekableIterator<VariantContext>(iterator2);
        this.comparator = new VariantContextComparator(dict);
    }

    @Override
    public boolean hasNext() {
        return this.iterator1.hasNext() || this.iterator2.hasNext();
    }

    @Override
    public VcTuple next() {
        if (!hasNext()) throw new IllegalStateException("next() called while hasNext() is false.");

        final VariantContext vc1 = this.iterator1.hasNext() ? this.iterator1.peek() : null;
        final VariantContext vc2 = this.iterator2.hasNext() ? this.iterator2.peek() : null;

        // If one or the other is null because there is no next, just return a one-sided tuple
        if (vc1 == null)       return new VcTuple(null,                  this.iterator2.next());
        else if (vc2 == null)  return new VcTuple(this.iterator1.next(), null);

        // Otherwise check the ordering and do the right thing
        final int ordering = this.comparator.compare(vc1, vc2);
        if (ordering == 0)     return new VcTuple(this.iterator1.next(), this.iterator2.next());
        else if (ordering < 0) return new VcTuple(this.iterator1.next(), null);
        else                   return new VcTuple(null,                  this.iterator2.next());
    }

    @Override public void remove() {
        throw new UnsupportedOperationException();
    }
}

class ConcordanceResults {
    private static class TwoState implements Comparable<TwoState> {
        final State state1, state2;

        private TwoState(final State state1, final State state2) {
            this.state1 = state1;
            this.state2 = state2;
        }

        @Override public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            return compareTo((TwoState) o) == 0;
        }

        @Override public int hashCode() {
            int result = state1.hashCode();
            result = 31 * result + state2.hashCode();
            return result;
        }

        @Override public int compareTo(final TwoState that) {
            int result = this.state1.compareTo(that.state1);
            if (result == 0) result = this.state2.compareTo(that.state2);
            return result;
        }
    }

    private final Histogram<TwoState> counts = new Histogram<TwoState>();
    private String sample1, sample2;

    ConcordanceResults(final String sample1, final String sample2) {
        this.sample1 = sample1;
        this.sample2 = sample2;
    }

    void add(final State state1, final State state2) {
        this.counts.increment(new TwoState(state1, state2));
    }

    long getCount(final State state1, final State state2) {
        final Histogram<TwoState>.Bin bin = this.counts.get(new TwoState(state1, state2));
        if (bin == null) return 0;
        else return (long) bin.getValue();
    }

    public String getSample2() { return sample2; }
    public String getSample1() { return sample1; }

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

    /** Sums the counts where the first state is contains in lhs and the second state is contained in rhs. */
    private double sum(final EnumSet<State> lhs, final EnumSet<State> rhs) {
        double result = 0;
        for (final State s1 : lhs) {
            for (final State s2 : rhs) {
                result += getCount(s1, s2);
            }
        }

        return result;
    }

}

