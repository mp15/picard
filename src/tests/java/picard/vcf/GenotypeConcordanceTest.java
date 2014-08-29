package picard.vcf;

import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

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



    @AfterClass
    public void teardown() {
        IOUtil.deleteDirectoryTree(OUTPUT_DATA_PATH);
    }

    @DataProvider(name = "genotypeConcordanceTestData")
    public Object[][] getGenotypeConcordanceTestData() {
        return new Object[][] {
                {CEU_TRIOS_SNPS_VCF, "NA12878", CEU_TRIOS_SNPS_VCF, "NA12878", new File(TEST_DATA_PATH, CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_GC)},
                {CEU_TRIOS_INDELS_VCF, "NA12878", CEU_TRIOS_INDELS_VCF, "NA12878", new File(TEST_DATA_PATH, CEU_TRIOS_INDELS_VS_CEU_TRIOS_INDELS_GC)},
                {CEU_TRIOS_SNPS_VCF, "NA12878", CEU_TRIOS_SNPS_FIRST_LINE_DIFF_VCF, "NA12878", new File(TEST_DATA_PATH, CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_FIRST_LINE_DIFF_GC)},
                {CEU_TRIOS_SNPS_VCF, "NA12878", CEU_TRIOS_SNPS_LAST_LINE_DIFF_VCF, "NA12878", new File(TEST_DATA_PATH, CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_LAST_LINE_DIFF_GC)},
                {CEU_TRIOS_SNPS_VCF, "NA12878", CEU_TRIOS_SNPS_DEL_LINE_VCF, "NA12878", new File(TEST_DATA_PATH, CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_DEL_LINE_GC)},
        };
    }

    @Test(dataProvider = "genotypeConcordanceTestData")
    public void testGenotypeConcordance(final File vcf1, final String sample1, final File vcf2, final String sample2, final File expectedOutputFile) throws Exception {
        final File outputFile = new File(OUTPUT_DATA_PATH, "actualGtConc.txt");
        outputFile.deleteOnExit();

        final GenotypeConcordance genotypeConcordance = new GenotypeConcordance();
        genotypeConcordance.VCF1 = vcf1;
        genotypeConcordance.SAMPLE1 = sample1;
        genotypeConcordance.VCF2 = vcf2;
        genotypeConcordance.SAMPLE2 = sample2;
        genotypeConcordance.OUTPUT = outputFile;
        final int returnCode = genotypeConcordance.instanceMain(new String[0]);
        Assert.assertEquals(returnCode, 0);

        IOUtil.assertFilesEqual(outputFile, expectedOutputFile);
    }


}
