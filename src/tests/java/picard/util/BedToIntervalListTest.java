package picard.util;

import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;

/**
 * @author nhomer
 */
public class BedToIntervalListTest {

    private final String TEST_DATA_DIR = "testdata/picard/util/BedToIntervalListTest";

    private String getFileContentsAsString(final File inputFile) throws IOException {
        final InputStream inputStream = IOUtil.openFileForReading(inputFile);
        final String inputBedString = IOUtil.readFully(inputStream);
        inputStream.close();
        return inputBedString;
    }

    private void doTest(final String inputBed, final String header) throws IOException {
        final File outputFile  = File.createTempFile("bed_to_interval_list_test.", ".interval_list");
        final BedToIntervalList program = new BedToIntervalList();
        final File inputBedFile = new File(TEST_DATA_DIR, inputBed);
        program.INPUT = inputBedFile;
        program.SEQUENCE_DICTIONARY = new File(TEST_DATA_DIR, header);
        program.OUTPUT = outputFile;
        program.doWork();

        // Assert they are equal
        Assert.assertEquals(getFileContentsAsString(new File(inputBedFile.getAbsolutePath() + ".interval_list")),
                getFileContentsAsString(outputFile));

        // remove the output file
        Assert.assertEquals(outputFile.delete(), true);
    }

    @Test(dataProvider = "testBedToIntervalListDataProvider")
    public void testBedToIntervalList(final String inputBed) throws IOException {
        doTest(inputBed, "header.sam");
    }

    @Test(dataProvider = "testBedToIntervalListOutOfBoundsDataProvider", expectedExceptions = PicardException.class)
    public void testBedToIntervalListOutOfBounds(final String inputBed) throws IOException {
        doTest(inputBed, "header.sam");
    }

    @DataProvider
    public Object[][] testBedToIntervalListDataProvider() {
        return new Object[][]{
                {"simple.bed"},
                {"overlapping.bed"}
        };
    }

    @DataProvider
    public Object[][] testBedToIntervalListOutOfBoundsDataProvider() {
        return new Object[][]{
                {"end_after_chr.bed"},
                {"end_before_chr.bed"},
                {"missing_chr.bed"},
                {"start_after_chr.bed"},
                {"start_before_chr.bed"}
        };
    }
}
