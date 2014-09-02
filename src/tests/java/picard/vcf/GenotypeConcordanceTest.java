package picard.vcf;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.FormatUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

public class GenotypeConcordanceTest {

    private static final File OUTPUT_DATA_PATH = IOUtil.createTempDir("GenotypeConcordanceTest", null);
    private static final File TEST_DATA_PATH = new File("testdata/picard/vcf/");

    // Test VCFs
    private static final File CEU_TRIOS_SNPS_VCF = new File(TEST_DATA_PATH, "CEUTrio-snps.vcf");
    private static final File CEU_TRIOS_INDELS_VCF = new File(TEST_DATA_PATH, "CEUTrio-indels.vcf");

    // Test that we notice a difference on the first line
    private static final File CEU_TRIOS_SNPS_FIRST_LINE_DIFF_VCF = new File(TEST_DATA_PATH, "CEUTrio-snps_first_line_diff.vcf");

    // Test that we notice a difference on the last line
    private static final File CEU_TRIOS_SNPS_LAST_LINE_DIFF_VCF = new File(TEST_DATA_PATH, "CEUTrio-snps_last_line_diff.vcf");

    // Test that we notice a deleted line
    private static final File CEU_TRIOS_SNPS_DEL_LINE_VCF = new File(TEST_DATA_PATH, "CEUTrio-snps_del_line.vcf");

    // Existing/expected results files
    private static final String CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_GC = "CEUTrio-snps_vs_CEUTrio-snps_GtConcordanceDiff.txt";
    private static final String CEU_TRIOS_INDELS_VS_CEU_TRIOS_INDELS_GC = "CEUTrio-indels_vs_CEUTrio-indels_GtConcordanceDiff.txt";
    private static final String CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_FIRST_LINE_DIFF_GC = "CEUTrio-snps_CEUTrio-snps_first_line_GtConcordanceDiff.txt";
    private static final String CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_LAST_LINE_DIFF_GC = "CEUTrio-snps_CEUTrio-snps_last_line_GtConcordanceDiff.txt";
    private static final String CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_DEL_LINE_GC = "CEUTrio-snps_CEUTrio-snps_del_line_GtConcordanceDiff.txt";
    private static final String CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_GC_ALL_ROWS = "CEUTrio-snps_vs_CEUTrio-snps_GtConcordanceDiff_AllRows.txt";
    private static final String CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_GC_MIN_GQ = "CEUTrio-snps_vs_CEUTrio-snps_GtConcordanceDiff_MinGq.txt";
    private static final String CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_GC_MIN_DP = "CEUTrio-snps_vs_CEUTrio-snps_GtConcordanceDiff_MinDp.txt";

    private static final File INTERVALS_FILE = new File(TEST_DATA_PATH, "IntervalListChr1Small.interval_list");

    @AfterClass
    public void teardown() {
        IOUtil.deleteDirectoryTree(OUTPUT_DATA_PATH);
    }

    @DataProvider(name = "genotypeConcordanceTestFileData")
    public Object[][] getGenotypeConcordanceTestFileData() {
        return new Object[][] {
                {CEU_TRIOS_SNPS_VCF, "NA12878", CEU_TRIOS_SNPS_VCF, "NA12878", null, null, false, new File(TEST_DATA_PATH, CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_GC)},
                {CEU_TRIOS_INDELS_VCF, "NA12878", CEU_TRIOS_INDELS_VCF, "NA12878", null, null, false, new File(TEST_DATA_PATH, CEU_TRIOS_INDELS_VS_CEU_TRIOS_INDELS_GC)},
                {CEU_TRIOS_SNPS_VCF, "NA12878", CEU_TRIOS_SNPS_FIRST_LINE_DIFF_VCF, "NA12878", null, null, false, new File(TEST_DATA_PATH, CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_FIRST_LINE_DIFF_GC)},
                {CEU_TRIOS_SNPS_VCF, "NA12878", CEU_TRIOS_SNPS_LAST_LINE_DIFF_VCF, "NA12878", null, null, false, new File(TEST_DATA_PATH, CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_LAST_LINE_DIFF_GC)},
                {CEU_TRIOS_SNPS_VCF, "NA12878", CEU_TRIOS_SNPS_DEL_LINE_VCF, "NA12878", null, null, false, new File(TEST_DATA_PATH, CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_DEL_LINE_GC)},
                {CEU_TRIOS_SNPS_VCF, "NA12878", CEU_TRIOS_SNPS_VCF, "NA12878", null, null, true, new File(TEST_DATA_PATH, CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_GC_ALL_ROWS)},
                {CEU_TRIOS_SNPS_VCF, "NA12878", CEU_TRIOS_SNPS_VCF, "NA12891", 40, null, false, new File(TEST_DATA_PATH, CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_GC_MIN_GQ)},
                {CEU_TRIOS_SNPS_VCF, "NA12878", CEU_TRIOS_SNPS_VCF, "NA12891", null, 40, false, new File(TEST_DATA_PATH, CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_GC_MIN_DP)},

        };
    }

    @Test(dataProvider = "genotypeConcordanceTestFileData")
    public void testGenotypeConcordance(final File vcf1, final String sample1, final File vcf2, final String sample2,
                                        final Integer minGq, final Integer minDp, final boolean outputAllRows,
                                        final File expectedOutputFile) throws Exception {
        final File outputFile = new File(OUTPUT_DATA_PATH, "actualGtConc.txt");
        outputFile.deleteOnExit();

        final GenotypeConcordance genotypeConcordance = new GenotypeConcordance();
        genotypeConcordance.TRUTH_VCF = vcf1;
        genotypeConcordance.TRUTH_SAMPLE = sample1;
        genotypeConcordance.CALL_VCF = vcf2;
        genotypeConcordance.CALL_SAMPLE = sample2;
        if (minGq != null)
            genotypeConcordance.MIN_GQ = minGq;
        if (minDp != null)
            genotypeConcordance.MIN_DP = minDp;
        genotypeConcordance.OUTPUT_ALL_ROWS = outputAllRows;
        genotypeConcordance.OUTPUT = outputFile;
        final int returnCode = genotypeConcordance.instanceMain(new String[0]);

        genotypeConcordance.getSnpCounter();
        Assert.assertEquals(returnCode, 0);

        IOUtil.assertFilesEqual(outputFile, expectedOutputFile);
    }

    @Test
    public void testGenotypeConcordanceDetails() throws Exception {
        final GenotypeConcordance genotypeConcordance = new GenotypeConcordance();
        genotypeConcordance.TRUTH_VCF = CEU_TRIOS_SNPS_VCF;
        genotypeConcordance.TRUTH_SAMPLE = "NA12878";
        genotypeConcordance.CALL_VCF = CEU_TRIOS_SNPS_VCF;
        genotypeConcordance.CALL_SAMPLE = "NA12878";
        int returnCode = genotypeConcordance.instanceMain(new String[0]);
        Assert.assertEquals(returnCode, 0);

        final Map<TwoState, Integer> nonZeroCounts = new HashMap<TwoState, Integer>();
        nonZeroCounts.put(new TwoState(State.HomRef, State.HomRef), 61);
        nonZeroCounts.put(new TwoState(State.Het, State.Het), 104);
        nonZeroCounts.put(new TwoState(State.HomVar, State.HomVar), 59);
        nonZeroCounts.put(new TwoState(State.FilteredVariant, State.FilteredVariant), 51);

        ConcordanceResults concordanceResults = genotypeConcordance.getSnpCounter();
        for (final State state1 : State.values()) {
            for (final State state2 : State.values()) {
                Integer expectedCount = nonZeroCounts.get(new TwoState(state1, state2));
                if (expectedCount == null) expectedCount = 0;
                Assert.assertEquals(concordanceResults.getCount(state1, state2), expectedCount.intValue());
            }
        }

        final FormatUtil fmt = new FormatUtil();

        Assert.assertEquals(fmt.format(concordanceResults.hetSensitivity()), "1");
        Assert.assertEquals(fmt.format(concordanceResults.hetPpv()), "1");
        Assert.assertEquals(fmt.format(concordanceResults.homVarSensitivity()), "1");
        Assert.assertEquals(fmt.format(concordanceResults.homVarPpv()), "1");
        Assert.assertEquals(fmt.format(concordanceResults.varSensitivity()), "1");
        Assert.assertEquals(fmt.format(concordanceResults.varPpv()), "1");

        // Now run it again with different samples
        genotypeConcordance.TRUTH_VCF = CEU_TRIOS_SNPS_VCF;
        genotypeConcordance.TRUTH_SAMPLE = "NA12878";
        genotypeConcordance.CALL_VCF = CEU_TRIOS_SNPS_VCF;
        genotypeConcordance.CALL_SAMPLE = "NA12891";
        returnCode = genotypeConcordance.instanceMain(new String[0]);
        Assert.assertEquals(returnCode, 0);

        nonZeroCounts.clear();
        nonZeroCounts.put(new TwoState(State.HomRef, State.HomRef), 30);
        nonZeroCounts.put(new TwoState(State.HomRef, State.Het), 31);
        nonZeroCounts.put(new TwoState(State.Het, State.HomRef), 30);
        nonZeroCounts.put(new TwoState(State.Het, State.Het), 50);
        nonZeroCounts.put(new TwoState(State.Het, State.HomVar), 24);
        nonZeroCounts.put(new TwoState(State.HomVar, State.Het), 18);
        nonZeroCounts.put(new TwoState(State.HomVar, State.HomVar), 41);
        nonZeroCounts.put(new TwoState(State.FilteredVariant, State.FilteredVariant), 51);

        concordanceResults = genotypeConcordance.getSnpCounter();
        for (final State state1 : State.values()) {
            for (final State state2 : State.values()) {
                Integer expectedCount = nonZeroCounts.get(new TwoState(state1, state2));
                if (expectedCount == null) expectedCount = 0;
                Assert.assertEquals(concordanceResults.getCount(state1, state2), expectedCount.intValue());
            }
        }

        Assert.assertEquals(fmt.format(concordanceResults.hetSensitivity()), "0.480769");
        Assert.assertEquals(fmt.format(concordanceResults.hetPpv()), "0.505051");
        Assert.assertEquals(fmt.format(concordanceResults.homVarSensitivity()), "0.694915");
        Assert.assertEquals(fmt.format(concordanceResults.homVarPpv()), "0.630769");
        Assert.assertEquals(fmt.format(concordanceResults.varSensitivity()), "0.815951");
        Assert.assertEquals(fmt.format(concordanceResults.varPpv()), "0.810976");
    }

    @Test
    public void testGenotypeConcordanceDetailsWithIntervals() throws Exception {
        final GenotypeConcordance genotypeConcordance = new GenotypeConcordance();
        genotypeConcordance.TRUTH_VCF = CEU_TRIOS_SNPS_VCF;
        genotypeConcordance.TRUTH_SAMPLE = "NA12878";
        genotypeConcordance.CALL_VCF = CEU_TRIOS_SNPS_VCF;
        genotypeConcordance.CALL_SAMPLE = "NA12878";
        genotypeConcordance.INTERVALS = Collections.singletonList(INTERVALS_FILE);

        int returnCode = genotypeConcordance.instanceMain(new String[0]);
        Assert.assertEquals(returnCode, 0);

        final Map<TwoState, Integer> nonZeroCounts = new HashMap<TwoState, Integer>();
        nonZeroCounts.put(new TwoState(State.HomRef, State.HomRef), 2);
        nonZeroCounts.put(new TwoState(State.Het, State.Het), 1);
        nonZeroCounts.put(new TwoState(State.FilteredVariant, State.FilteredVariant), 2);

        ConcordanceResults concordanceResults = genotypeConcordance.getSnpCounter();
        for (final State state1 : State.values()) {
            for (final State state2 : State.values()) {
                Integer expectedCount = nonZeroCounts.get(new TwoState(state1, state2));
                if (expectedCount == null) expectedCount = 0;
                Assert.assertEquals(concordanceResults.getCount(state1, state2), expectedCount.intValue());
            }
        }

        final FormatUtil fmt = new FormatUtil();

        Assert.assertEquals(fmt.format(concordanceResults.hetSensitivity()), "1");
        Assert.assertEquals(fmt.format(concordanceResults.hetPpv()), "1");
        Assert.assertEquals(fmt.format(concordanceResults.homVarSensitivity()), "?");
        Assert.assertEquals(fmt.format(concordanceResults.homVarPpv()), "?");
        Assert.assertEquals(fmt.format(concordanceResults.varSensitivity()), "1");
        Assert.assertEquals(fmt.format(concordanceResults.varPpv()), "1");

        // Now run it again with different samples
        genotypeConcordance.TRUTH_VCF = CEU_TRIOS_SNPS_VCF;
        genotypeConcordance.TRUTH_SAMPLE = "NA12878";
        genotypeConcordance.CALL_VCF = CEU_TRIOS_SNPS_VCF;
        genotypeConcordance.CALL_SAMPLE = "NA12891";
        genotypeConcordance.INTERVALS = Collections.singletonList(INTERVALS_FILE);
        returnCode = genotypeConcordance.instanceMain(new String[0]);
        Assert.assertEquals(returnCode, 0);

        nonZeroCounts.clear();
        nonZeroCounts.put(new TwoState(State.HomRef, State.HomRef), 1);
        nonZeroCounts.put(new TwoState(State.HomRef, State.Het), 1);
        nonZeroCounts.put(new TwoState(State.Het, State.Het), 1);
        nonZeroCounts.put(new TwoState(State.FilteredVariant, State.FilteredVariant), 2);

        concordanceResults = genotypeConcordance.getSnpCounter();
        for (final State state1 : State.values()) {
            for (final State state2 : State.values()) {
                Integer expectedCount = nonZeroCounts.get(new TwoState(state1, state2));
                if (expectedCount == null) expectedCount = 0;
                Assert.assertEquals(concordanceResults.getCount(state1, state2), expectedCount.intValue());
            }
        }

        Assert.assertEquals(fmt.format(concordanceResults.hetSensitivity()), "1");
        Assert.assertEquals(fmt.format(concordanceResults.hetPpv()), "0.5");
        Assert.assertEquals(fmt.format(concordanceResults.homVarSensitivity()), "?");
        Assert.assertEquals(fmt.format(concordanceResults.homVarPpv()), "?");
        Assert.assertEquals(fmt.format(concordanceResults.varSensitivity()), "1");
        Assert.assertEquals(fmt.format(concordanceResults.varPpv()), "0.5");
    }

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
}
