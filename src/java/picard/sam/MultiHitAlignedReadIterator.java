/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package picard.sam;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.FastBAMRecord;
import htsjdk.samtools.ReadRecord;
import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.filter.FilteringIterator;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.PeekableIterator;
import picard.PicardException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.NoSuchElementException;


/**
 * Iterate over queryname-sorted SAM, and return each group of reads with the same queryname.  Unmapped reads
 * are filtered out, as are alignments that don't seem to match any part of the reference.
 * If there are multiple hits for the same read, and the first and second ends need to be correlated,
 * then they are sorted by hit index. Supplemental alignments are discarded, with a logged message.
 * A set of hits for a single query may then be filtered with a caller-supplied filter, which will remove any
 * alignments that do not pass the filter.  If the primary alignment is removed, the best-mapping secondary alignment
 * or alignment pair will be marked as primary.
 *
 *
 * @throws IllegalStateException if the input is not queryname-sorted.
 */
class MultiHitAlignedReadIterator implements CloseableIterator<HitsForInsert> {
    private final PeekableIterator<ReadRecord> peekIterator;
    private final SAMRecordQueryNameComparator queryNameComparator = new SAMRecordQueryNameComparator();
    private final PrimaryAlignmentSelectionStrategy primaryAlignmentSelectionStrategy;

    private HitsForInsert theNext = null;

    /**
     *
     * @param querynameOrderIterator
     * @param primaryAlignmentSelectionStrategy Algorithm for selecting primary alignment when it is not clear from
     *                                          the input what should be primary.
     */
    MultiHitAlignedReadIterator(final CloseableIterator<ReadRecord> querynameOrderIterator,
                                final PrimaryAlignmentSelectionStrategy primaryAlignmentSelectionStrategy) {
        this.primaryAlignmentSelectionStrategy = primaryAlignmentSelectionStrategy;
        peekIterator = new PeekableIterator<ReadRecord>(new FilteringIterator(querynameOrderIterator,
                new SamRecordFilter() {
                    // Filter unmapped reads.
                    public boolean filterOut(final ReadRecord record) {
                        return record.getReadUnmappedFlag() || SAMUtils.cigarMapsNoBasesToRef(record);
                    }
                    public boolean filterOut(final ReadRecord first, final ReadRecord second) {
                        return ((first.getReadUnmappedFlag() || SAMUtils.cigarMapsNoBasesToRef(first))
                                && (second.getReadUnmappedFlag() || SAMUtils.cigarMapsNoBasesToRef(second)));
                    }
                }));


        advance();
    }

    public void close() {
        peekIterator.close();
    }

    public boolean hasNext() {
        return theNext != null;
    }

    /**
     * @throws IllegalStateException if the input is not queryname-sorted.
     */
    public HitsForInsert next() {
        if (!hasNext()) throw new NoSuchElementException();
        final HitsForInsert ret = theNext;
        advance();
        return ret;
    }

    private void advance() {
        while (peekIterator.hasNext()) {
            theNext = nextMaybeEmpty();
            if (theNext.numHits() > 0) return;
        }
        theNext = null;
    }

    private HitsForInsert nextMaybeEmpty() {
        if (!peekIterator.hasNext()) throw new IllegalStateException();
        final String readName = peekIterator.peek().getReadName();
        final HitsForInsert hits = new HitsForInsert();

        Boolean isPaired = null;

        // Accumulate the alignments matching readName.
        do {
            final ReadRecord rec = peekIterator.next();
            replaceHardWithSoftClips((FastBAMRecord) rec);
            // It is critical to do this here, because SamAlignmentMerger uses this exception to determine
            // if the aligned input needs to be sorted.
            if (peekIterator.hasNext() && queryNameComparator.fileOrderCompare(rec, peekIterator.peek()) > 0) {
                throw new IllegalStateException("Underlying iterator is not queryname sorted: " +
                rec + " > " + peekIterator.peek());
            }

            if (isPaired == null) {
                isPaired = rec.getReadPairedFlag();
            } else if (isPaired != rec.getReadPairedFlag()) {
                throw new PicardException("Got a mix of paired and unpaired alignments for read " + readName);
            }

            // Records w/ a supplemental flag are stashed to the side until the primary alignment has
            // been determined, and then re-added into the process later
            if (!rec.getReadPairedFlag() || rec.getFirstOfPairFlag()) {
                if (rec.getSupplementaryAlignmentFlag()) {
                    hits.addSupplementalFirstOfPairOrFragment(rec);
                } else {
                    hits.addFirstOfPairOrFragment(rec);
                }
            } else if (rec.getSecondOfPairFlag()) {
                if (rec.getSupplementaryAlignmentFlag()) {
                    hits.addSupplementalSecondOfPair(rec);
                } else {
                    hits.addSecondOfPair(rec);
                }
            } else throw new PicardException("Read is marked as pair but neither first or second: " + readName);
        } while (peekIterator.hasNext() && peekIterator.peek().getReadName().equals(readName));

        // If there is no more than one alignment for each end, no need to do any coordination.
        if (hits.numHits() <= 1) {
            // No HI tags needed if only a single hit
            if (hits.getFirstOfPair(0) != null) {
                hits.getFirstOfPair(0).setAttribute(SAMTag.HI.name(), null);
                hits.getFirstOfPair(0).setNotPrimaryAlignmentFlag(false);
            }
            if (hits.getSecondOfPair(0) != null) {
                hits.getSecondOfPair(0).setAttribute(SAMTag.HI.name(), null);
                hits.getSecondOfPair(0).setNotPrimaryAlignmentFlag(false);
            }
        } else {
            primaryAlignmentSelectionStrategy.pickPrimaryAlignment(hits);
        }

        // Used to check that alignments for first and second were correlated, but this is no longer required.
        return hits;
    }

    /** Replaces hard clips with soft clips and fills in bases and qualities with dummy values as needed. */
    private void replaceHardWithSoftClips(final FastBAMRecord rec) {
        if (rec.getReadUnmappedFlag()) return;
        if (rec.getCigarLength() == 0) return;

        /*
        List<CigarElement> elements = rec.getCigar().getCigarElements();
        final CigarElement first = elements.get(0);
        final CigarElement last  = elements.size() == 1 ? null : elements.get(elements.size()-1);
        final int startHardClip = first.getOperator() == CigarOperator.H ? first.getLength() : 0;
        final int endHardClip   = (last != null && last.getOperator() == CigarOperator.H) ? last.getLength() : 0;

         */

        final CigarOperator firstOp = rec.getCigarOp(0);
        final CigarOperator lastOp  = rec.getCigarLength() == 1 ? null : rec.getCigarOp(rec.getCigarLength()-1);
        final int startHardClip = firstOp == CigarOperator.H ? rec.getCigarOpLength(0) : 0;
        final int endHardClip   = (lastOp != null && lastOp == CigarOperator.H) ? rec.getCigarOpLength(rec.getCigarLength()-1) : 0;

        if (startHardClip + endHardClip > 0) {
            final int len = rec.getReadBases().length + startHardClip + endHardClip;

            // Fix the basecalls
            final byte[] bases = new byte[len];
            Arrays.fill(bases, (byte) 'N');
            System.arraycopy(rec.getReadBases(), 0, bases, startHardClip, rec.getReadBases().length);

            // Fix the quality scores
            final byte[] quals = new byte[len];
            Arrays.fill(quals, (byte) 2  );
            System.arraycopy(rec.getBaseQualities(), 0, quals, startHardClip, rec.getBaseQualities().length);

            // Fix the cigar!
            if (startHardClip > 0) rec.setCigarOp(0, CigarOperator.S);
            if (endHardClip   > 0) rec.setCigarOp(rec.getCigarLength()-1, CigarOperator.S);

            // Set the update structures on the new record
            rec.setReadBases(bases);
            rec.setBaseQualities(quals);
        }
   }

    /** Unsupported operation. */
    public void remove() {
        throw new UnsupportedOperationException();
    }
}
