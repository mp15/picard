/*
 * The MIT License
 *
 * Copyright (c) 2012 The Broad Institute
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
package picard.sam.markduplicates;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.Closeable;
import java.io.File;
import java.util.ArrayList;

public abstract class AbstractMarkDuplicateFindingAlgorithmTest {

    protected abstract AbstractMarkDuplicateFindingAlgorithmTester getTester();

    protected final static int DEFAULT_BASE_QUALITY = 10;

    @Test
    public void testSingleUnmappedFragment() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addUnmappedFragment(-1, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testSingleUnmappedPair() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addUnmappedPair(-1, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testSingleMappedFragment() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedFragment(1, 1, false, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testTwoMappedFragments() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedFragment(0, 1, false, DEFAULT_BASE_QUALITY);
        tester.addMappedFragment(0, 1, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.runTest();
    }

    @Test
    public void testSingleMappedPair() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(1, 1, 100, false, false, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testSingleMappedFragmentAndSingleMappedPair() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedFragment(1, 1, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(1, 1, 100, false, false, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairs() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(1, 1, 100, false, false, DEFAULT_BASE_QUALITY);
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.runTest();
    }

    @Test
    public void testThreeMappedPairs() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(1, 1, 100, false, false, DEFAULT_BASE_QUALITY);
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.runTest();
    }

    @Test
    public void testSingleMappedFragmentAndTwoMappedPairs() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedFragment(1, 1, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(1, 1, 100, false, false, DEFAULT_BASE_QUALITY);
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairsAndTerminalUnmappedFragment() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(1, 1, 100, false, false, DEFAULT_BASE_QUALITY);
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addUnmappedFragment(-1, DEFAULT_BASE_QUALITY); // unmapped fragment at end of file
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairsAndTerminalUnmappedPair() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(1, 1, 100, false, false, DEFAULT_BASE_QUALITY);
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addUnmappedPair(-1, DEFAULT_BASE_QUALITY); // unmapped pair at end of file
        tester.runTest();
    }

    @Test
    public void testOpticalDuplicateFinding() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();

        // explicitly creating 1 expected optical duplicate pair
        tester.setExpectedOpticalDuplicate(1);

        // pass in the read names manually, in order to control duplicates vs optical duplicates
        tester.addMatePair("READ0:1:1:1:1", 1, 1, 100, false, false, false, false, "50M", "50M", false, true, false,
                           false, false, DEFAULT_BASE_QUALITY); // non-duplicate mapped pair to start
        tester.addMatePair("READ1:1:1:1:300", 1, 1, 100, false, false, true, true, "50M", "50M", false, true, false,
                           false, false, DEFAULT_BASE_QUALITY); // duplicate pair, NOT optical duplicate (delta-Y > 100)
        tester.addMatePair("READ2:1:1:1:50", 1, 1, 100, false, false, true, true, "50M", "50M", false, true, false,
                           false, false, DEFAULT_BASE_QUALITY); // duplicate pair, expected optical duplicate (delta-X and delta-Y < 100)
        tester.runTest();
    }

    @Test
    public void testOpticalDuplicateClusterSamePosition() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.setExpectedOpticalDuplicate(0);
        tester.addMatePair("RUNID:7:1203:2886:82292", 1, 485253, 485253, false, false, true, true, "42M59S", "59S42M", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.addMatePair("RUNID:7:1203:2884:16834", 1, 485253, 485253, false, false, false, false, "59S42M", "42M59S", true, false, false, false, false, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testOpticalDuplicateClusterOneEndSamePositionNoCluster() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.setExpectedOpticalDuplicate(1);
        tester.addMatePair("RUNID:7:2205:17939:39728", 1, 485328, 485312, false, false, false, false, "55M46S", "30S71M", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.addMatePair("RUNID:7:2205:17949:39745", 1, 485328, 485328, false, false, true, true, "55M46S", "46S55M", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairsAndMappedSecondaryFragment() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(1, 1, 100, false, false, DEFAULT_BASE_QUALITY);
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedFragment(1, 200, false, DEFAULT_BASE_QUALITY, true); // mapped non-primary fragment
        tester.runTest();
    }

    @Test
    public void testMappedFragmentAndMappedPairFirstOfPairNonPrimary() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedFragment(1, 1, false,DEFAULT_BASE_QUALITY);
        tester.addMatePair(1, 200, 0, false, true, false, false, "54M22S", null, false, false, true, true, false, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairsMatesSoftClipped() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(1, 10022, 10051, false, false, "76M", "8S68M", false, true, false, DEFAULT_BASE_QUALITY);
        tester.addMappedPair(1, 10022, 10063, false, false, "76M", "5S71M", false, true, false, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairsWithSoftClipping() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        // NB: no duplicates
        // 5'1: 2, 5'2:46+73M=118
        // 5'1: 2, 5'2:51+68M=118
        tester.addMappedPair(1, 2, 46, false, false, "6S42M28S", "3S73M", false, DEFAULT_BASE_QUALITY);
        tester.addMappedPair(1, 2, 51, true, true, "6S42M28S", "8S68M", false, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairsWithSoftClippingFirstOfPairOnlyNoMateCigar() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.setNoMateCigars(true);
        // NB: no duplicates
        // 5'1: 2, 5'2:46+73M=118
        // 5'1: 2, 5'2:51+68M=118
        tester.addMappedPair(1, 12, 46, false, false, "6S42M28S", null, true, DEFAULT_BASE_QUALITY); // only add the first one
        tester.addMappedPair(1, 12, 51, false, false, "6S42M28S", null, true, DEFAULT_BASE_QUALITY); // only add the first one
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairsWithSoftClippingBoth() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(1, 10046, 10002, false, false, "3S73M", "6S42M28S", true, false, false, DEFAULT_BASE_QUALITY);
        tester.addMappedPair(1, 10051, 10002, true, true, "8S68M", "6S48M22S", true, false, false, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testMatePairSecondUnmapped() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMatePair(1, 10049, 10049, false, true, false, false, "11M2I63M", null, false, false, false, false, false, DEFAULT_BASE_QUALITY);   // neither are duplicates
        tester.runTest();
    }

    @Test
    public void testMatePairFirstUnmapped() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMatePair(1, 10056, 10056, true, false, false, false, null, "54M22S", false, false, false, false, false, DEFAULT_BASE_QUALITY);    // neither are duplicates
        tester.runTest();
    }

    @Test
    public void testMappedFragmentAndMatePairSecondUnmapped() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMatePair(1, 10049, 10049, false, true, false, false, "11M2I63M", null, false, false, false, false, false, DEFAULT_BASE_QUALITY);
        tester.addMappedFragment(1, 10049, true, DEFAULT_BASE_QUALITY); // duplicate
        tester.runTest();
    }

    @Test
    public void testMappedFragmentAndMatePairFirstUnmapped() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMatePair(1, 10049, 10049, true, false, false, false, null, "11M2I63M", false, false, false, false, false, DEFAULT_BASE_QUALITY);
        tester.addMappedFragment(1, 10049, true, DEFAULT_BASE_QUALITY); // duplicate
        tester.runTest();
    }

    @Test
    public void testMappedPairAndMatePairSecondUnmapped() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMatePair(1, 10040, 10040, false, true, true, false, "76M", null, false, false, false, false, false, DEFAULT_BASE_QUALITY); // second a duplicate,
        // second end unmapped
        tester.addMappedPair(1, 10189, 10040, false, false, "41S35M", "65M11S", true, false, false, DEFAULT_BASE_QUALITY); // mapped OK
        tester.runTest();
    }

    @Test
    public void testMappedPairAndMatePairFirstUnmapped() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMatePair(1, 10040, 10040, true, false, false, true,  null, "76M", false, false, false, false, false, DEFAULT_BASE_QUALITY); // first a duplicate,
        // first end unmapped
        tester.addMappedPair(1, 10189, 10040, false, false, "41S35M", "65M11S", true, false, false, DEFAULT_BASE_QUALITY); // mapped OK
        tester.runTest();
    }

    // TODO: fails on MarkDuplicatesWithMateCigar
    @Test
    public void testMappedPairAndMatePairFirstOppositeStrandSecondUnmapped() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        // first end mapped OK -, second end unmapped
        tester.addMatePair(1, 484071, 484071, false, true, false, false,  "66S35M", null, true, false, false, false, false, DEFAULT_BASE_QUALITY);
        // mapped OK +/-
        tester.addMappedPair(1, 484105, 484075, false, false, "35M66S", "30S71M", false, true, false, DEFAULT_BASE_QUALITY);
        tester.runTest();
    }

    @Test
    public void testMappedPairAndMappedFragmentAndMatePairSecondUnmapped() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMatePair(1, 10040, 10040, false, true, true, false, "76M", null, false, false, false, false, false, DEFAULT_BASE_QUALITY); // first a duplicate,
        // second end unmapped
        tester.addMappedPair(1, 10189, 10040, false, false, "41S35M", "65M11S", true, false, false, DEFAULT_BASE_QUALITY); // mapped OK
        tester.addMappedFragment(1, 10040, true, DEFAULT_BASE_QUALITY); // duplicate
        tester.runTest();
    }

    @Test
    public void testMappedPairAndMappedFragmentAndMatePairFirstUnmapped() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMatePair(1, 10040, 10040, true, false, false, true, null, "76M", false, false, false, false, false, DEFAULT_BASE_QUALITY); // first a duplicate,
        // first end unmapped
        tester.addMappedPair(1, 10189, 10040, false, false, "41S35M", "65M11S", true, false, false, DEFAULT_BASE_QUALITY); // mapped OK
        tester.addMappedFragment(1, 10040, true, DEFAULT_BASE_QUALITY); // duplicate
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairsWithOppositeOrientations() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(1, 10182, 10038, false, false, "32S44M", "66M10S", true, false, false, DEFAULT_BASE_QUALITY); // -/+
        tester.addMappedPair(1, 10038, 10182, true, true, "70M6S", "32S44M", false, true, false, DEFAULT_BASE_QUALITY); // +/-, both are duplicates
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairsWithOppositeOrientationsNumberTwo() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(1, 10038, 10182, false, false, "70M6S", "32S44M", false, true, false, DEFAULT_BASE_QUALITY); // +/-, both are duplicates
        tester.addMappedPair(1, 10182, 10038, true, true, "32S44M", "66M10S", true, false, false, DEFAULT_BASE_QUALITY); // -/+
        tester.runTest();
    }

    @Test
    public void testThreeMappedPairsWithMatchingSecondMate() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        // Read0 and Read2 are duplicates
        // 10181+35=10216, 10058
        tester.addMappedPair(1, 10181, 10058, false, false, "41S35M", "47M29S", true, false, false, DEFAULT_BASE_QUALITY); // -/+
        // 10181+37=10218, 10058
        tester.addMappedPair(1, 10181, 10058, false, false, "37S39M", "44M32S", true, false, false, DEFAULT_BASE_QUALITY); // -/+
        // 10180+36=10216, 10058
        tester.addMappedPair(1, 10180, 10058, true, true, "36S40M", "50M26S", true, false, false, DEFAULT_BASE_QUALITY); // -/+, both are duplicates
        tester.runTest();
    }

    @Test
    public void testMappedPairWithSamePosition() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(1, 4914, 4914, false, false, "37M39S", "73M3S", false, false, false, DEFAULT_BASE_QUALITY); // +/+
        tester.runTest();
    }

    @Test
    public void testMappedPairWithSamePositionSameCigar() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(1, 4914, 4914, false, false, "37M39S", "37M39S", false, false, false, DEFAULT_BASE_QUALITY); // +/+
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairWithSamePosition() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(0, 5604914, 5604914, false, false, "37M39S", "73M3S", false, false, false, DEFAULT_BASE_QUALITY); // +/+
        tester.addMappedPair(0, 5604914, 5604914, true, true, "37M39S", "73M3S", false, false, false, DEFAULT_BASE_QUALITY); // +/+
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairWithSamePositionDifferentStrands() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(0, 5604914, 5604914, false, false, "50M", "50M", true, false, false, DEFAULT_BASE_QUALITY); // +/-
        tester.addMappedPair(0, 5604914, 5604914, true, true, "50M", "50M", false, true, false, DEFAULT_BASE_QUALITY); // -/+
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairWithSamePositionDifferentStrands2() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(0, 5604914, 5604915, false, false, "50M", "50M", true, false, false, DEFAULT_BASE_QUALITY); // +/-
        tester.addMappedPair(0, 5604915, 5604914, true, true, "50M", "50M", false, true, false, DEFAULT_BASE_QUALITY); // -/+
        tester.runTest();
    }

    @Test
    public void testMappedPairWithFirstEndSamePositionAndOther() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(0, 5604914, 5605914, false, false, "37M39S", "73M3S", false, false, false, DEFAULT_BASE_QUALITY); // +/+
        tester.addMappedPair(0, 5604914, 5604914, false, false, "37M39S", "73M3S", false, false, false, DEFAULT_BASE_QUALITY); // +/+
        tester.runTest();
    }

    @Test
    public void testTwoGroupsOnDifferentChromosomesOfTwoFragments() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedFragment(0, 1, false, DEFAULT_BASE_QUALITY);
        tester.addMappedFragment(0, 1, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedFragment(1, 1, false, DEFAULT_BASE_QUALITY);
        tester.addMappedFragment(1, 1, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.runTest();
    }


    @Test
    public void testTwoGroupsOnDifferentChromosomesOfTwoMappedPairs() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(0, 1, 100, false, false, DEFAULT_BASE_QUALITY);
        tester.addMappedPair(0, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(1, 1, 100, false, false, DEFAULT_BASE_QUALITY);
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.runTest();
    }

    @Test
    public void testTwoGroupsOnDifferentChromosomesOfThreeMappedPairs() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(0, 1, 100, false, false, DEFAULT_BASE_QUALITY);
        tester.addMappedPair(0, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(0, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(1, 1, 100, false, false, DEFAULT_BASE_QUALITY);
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.runTest();
    }

    @Test
    public void testThreeGroupsOnDifferentChromosomesOfThreeMappedPairs() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(0, 1, 100, false, false, DEFAULT_BASE_QUALITY);
        tester.addMappedPair(0, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(0, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(1, 1, 100, false, false, DEFAULT_BASE_QUALITY);
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(2, 1, 100, false, false , DEFAULT_BASE_QUALITY);
        tester.addMappedPair(2, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.addMappedPair(2, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.runTest();
    }

    @Test
    public void testBulkFragmentsNoDuplicates() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        for(int position = 1; position <= 10000; position += 1) {
            tester.addMappedFragment(0, position, false, "100M", DEFAULT_BASE_QUALITY);
        }
        tester.runTest();
    }

    @Test
    public void testBulkFragmentsWithDuplicates() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        for(int position = 1; position <= 10000; position += 1) {
            tester.addMappedFragment(0, position, false, "100M", DEFAULT_BASE_QUALITY);
            tester.addMappedFragment(0, position, true, "100M", DEFAULT_BASE_QUALITY);
            tester.addMappedFragment(0, position, true, "100M", DEFAULT_BASE_QUALITY);
            tester.addMappedFragment(0, position, true, "100M", DEFAULT_BASE_QUALITY);
            tester.addMappedFragment(0, position, true, "100M", DEFAULT_BASE_QUALITY);
        }
        tester.runTest();
    }

    @Test
    public void testStackOverFlowPairSetSwap() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();

        File input = new File("testdata/picard/sam/MarkDuplicates/markDuplicatesWithMateCigar.pairSet.swap.sam");
        SamReader reader = SamReaderFactory.makeDefault().open(input);
        tester.setHeader(reader.getFileHeader());
        for (final SAMRecord record : reader) {
            tester.addRecord(record);
        }
        CloserUtil.close(reader);
        tester.setExpectedOpticalDuplicate(1);
        tester.runTest();
    }

    @Test
    public void testSecondEndIsBeforeFirstInCoordinate() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.addMappedPair(0, 108855339, 108855323, false, false, "33S35M", "17S51M", true, true, false, DEFAULT_BASE_QUALITY); // +/-
        tester.runTest();
    }

    @Test
    public void testPathologicalOrderingAtTheSamePosition() {
        final AbstractMarkDuplicateFindingAlgorithmTester tester = getTester();
        tester.setExpectedOpticalDuplicate(1);
        tester.addMatePair("RUNID:3:1:15013:113051", 0, 129384554, 129384554, false, false, false, false, "68M", "68M", false, false, false, false, false, DEFAULT_BASE_QUALITY);
        tester.addMatePair("RUNID:3:1:15029:113060", 0, 129384554, 129384554, false, false, true, true, "68M", "68M", false, false, false, false, false, DEFAULT_BASE_QUALITY);

        // Create the pathology
        CloseableIterator<SAMRecord> iterator = tester.getRecordIterator();
        int[] qualityOffset = {20, 30, 10, 40}; // creates an interesting pathological ordering
        int index = 0;
        while (iterator.hasNext()) {
            final SAMRecord record = iterator.next();
            byte[] quals = new byte[record.getReadLength()];
            for (int i = 0; i < record.getReadLength(); i++) {
                quals[i] = (byte)(qualityOffset[index] + 10);
            }
            record.setBaseQualities(quals);
            index++;
        }
        iterator.close();

        // Run the test
        tester.runTest();
    }
}
