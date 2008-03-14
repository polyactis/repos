/*
 * FusionCELDataTest.java
 * JUnit based test
 *
 * Created on October 17, 2005, 4:55 PM
 */

package affymetrix.fusion.cel;

import java.io.File;
import java.util.Arrays;
import java.util.List;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import affymetrix.fusion.FusionTagValuePair;

/**
 * 
 * @author ljevon
 */
public class FusionCELDataTest extends TestCase {

	public FusionCELDataTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(FusionCELDataTest.class);

		return suite;
	}

	/**
	 * Test of getError method, of class affymetrix.gcos.fusion.FusionCELData.
	 */
	public void testGetError() {
	}

	/**
	 * Test of getFileName method, of class affymetrix.gcos.cel.FusionCELData.
	 */
	public void testFileName() {
		FusionCELData c = new FusionCELData();
		c.setFileName("asdf");
		assertEquals(c.getFileName(), "asdf");
	}

	public void testFileId_from_GCOS() throws Exception {
		FusionCELData c = new FusionCELData();
		c.setFileName(new File("..\\..\\..\\CPPTest\\data\\test.v4.CEL").getCanonicalPath());
		assertTrue(c.readHeader());
		assertEquals(c.getFileId(), null);
	}

	public void testFileId_from_Calvin() throws Exception {
		FusionCELData c = new FusionCELData();
		c.setFileName(new File("..\\..\\..\\..\\calvin_files\\fusion\\data\\small_cel_file_partial_datheader")
				.getCanonicalPath());
		assertTrue(c.readHeader());
		String s = "0000065535-1147280233-0000005844-0000011008-0000006224";
		assertTrue(Arrays.equals(c.getFileId().getGuid(), s.getBytes("US-ASCII")));
	}

	public void testLibPackage_from_Calvin() throws Exception {
		FusionCELData c = new FusionCELData();
		c.setFileName(new File("..\\..\\..\\..\\calvin_files\\fusion\\data\\small_cel_file_partial_datheader")
				.getCanonicalPath());
		assertTrue(c.readHeader());
		assertEquals(c.getLibraryPackageName().compareTo("Hg-small-lib-package"), 0);
	}

	public void testLibPackage_from_GCOS() throws Exception {
		FusionCELData c = new FusionCELData();
		c.setFileName(new File("..\\..\\..\\CPPTest\\data\\test.v4.CEL").getCanonicalPath());
		assertTrue(c.readHeader());
		assertEquals(c.getLibraryPackageName(), null);
	}

	public void testMasterFile_from_Calvin() throws Exception {
		FusionCELData c = new FusionCELData();
		c.setFileName(new File("..\\..\\..\\..\\calvin_files\\fusion\\data\\small_cel_file_partial_datheader")
				.getCanonicalPath());
		assertTrue(c.readHeader());
		assertEquals(c.getMasterFileName().compareTo("Hg-small-master-file"), 0);
	}

	public void testMasterFile_from_GCOS() throws Exception {
		FusionCELData c = new FusionCELData();
		c.setFileName(new File("..\\..\\..\\CPPTest\\data\\test.v4.CEL").getCanonicalPath());
		assertTrue(c.readHeader());
		assertEquals(c.getMasterFileName(), null);
	}

	private void checkHeader(FusionCELData c, boolean checkMaskOutliers) throws Exception {
		assertEquals(c.getRows(), 126);
		assertEquals(c.getCols(), 126);
		assertEquals(c.getCells(), 126 * 126);
		if (checkMaskOutliers) {
			assertEquals(c.getNumMasked(), 1);
			assertEquals(c.getNumOutliers(), 3917);
		}

		List<FusionTagValuePair> p = c.getParameters();
		String[] ptags = { "Percentile", "CellMargin", "OutlierHigh", "OutlierLow", "AlgVersion", "FixedCellSize",
				"IgnoreOutliersInShiftRows", "FeatureExtraction", "UseSubgrids", "RandomizePixels", "ErrorBasis", "StdMult" };
		String[] pvalues = { "75", "2", "1.500", "1.004", "6.0", "FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "StdvMean",
				"1.000000" };

		assertEquals(p.size(), pvalues.length);
		for (int i = 0; i < p.size(); i++) {
			FusionTagValuePair param = p.get(i);
			assertEquals(param.getTag(), ptags[i]);
			assertEquals(param.getValue(), pvalues[i]);
		}

		assertEquals(c.getAlg(), "Percentile");
		assertEquals(c.getChipType(), "Test3");
		assertEquals(c.getCellMargin(), 2);
		assertEquals(c.getGridCorners().getUpperLeft().getX(), 239.0f, 0.00001f);
		assertEquals(c.getGridCorners().getUpperLeft().getY(), 234.0f, 0.00001f);
		assertEquals(c.getGridCorners().getUpperRight().getX(), 4504.0f, 0.00001f);
		assertEquals(c.getGridCorners().getUpperRight().getY(), 232.0f, 0.00001f);
		assertEquals(c.getGridCorners().getLowerRight().getX(), 4505.0f, 0.00001f);
		assertEquals(c.getGridCorners().getLowerRight().getY(), 4497.0f, 0.00001f);
		assertEquals(c.getGridCorners().getLowerLeft().getX(), 240.0f, 0.00001f);
		assertEquals(c.getGridCorners().getLowerLeft().getY(), 4500.0f, 0.00001f);
		assertEquals(
				c.getDatHeader(),
				"[81..46222]  Test3:CLS=4733 RWS=4733 XIN=3  YIN=3  VE=17        2.0 02/25/ 0 10:35:25       GridVerify=None  Test3.1sq                  6");
	}

	private void checkBody(FusionCELData c) throws Exception {
		FusionCELFileEntryType entry = new FusionCELFileEntryType();
		c.getEntry(0, entry);
		assertEquals(entry.getIntensity(), 17047.0f, 0.00001f);
		assertEquals(entry.getStdv(), 10412.4f, 0.00001f);
		assertEquals(entry.getPixels(), 1024);
		assertEquals(c.getIntensity(0), 17047.0f, 0.00001f);
		assertEquals(c.getStdv(0), 10412.4f, 0.00001f);
		assertEquals(c.getPixels(0), 1024);

		c.getEntry(127, entry);
		assertEquals(entry.getIntensity(), 1094.0f, 0.00001f);
		assertEquals(entry.getStdv(), 1407.8f, 0.00001f);
		assertEquals(entry.getPixels(), 1024);
		assertEquals(c.getIntensity(127), 1094.0f, 0.00001f);
		assertEquals(c.getStdv(127), 1407.8f, 0.00001f);
		assertEquals(c.getPixels(127), 1024);

		assertEquals(c.indexToX(0), 0);
		assertEquals(c.indexToY(0), 0);
		assertEquals(c.indexToX(127), 1);
		assertEquals(c.indexToY(127), 1);

		assertEquals(c.xyToIndex(0, 0), 0);
		assertEquals(c.xyToIndex(1, 1), 127);
		assertEquals(c.xyToIndex(1, 0), 1);
		assertEquals(FusionCELData.xyToIndex(0, 0, 126, 126), 0);
		assertEquals(FusionCELData.xyToIndex(1, 1, 126, 126), 127);
		assertEquals(FusionCELData.xyToIndex(1, 0, 126, 126), 1);

		assertTrue(c.isMasked(0));
		assertEquals(c.isMasked(1), false);
		assertTrue(c.isMasked(0, 0));
		assertEquals(c.isMasked(1, 0), false);

		assertTrue(c.isOutlier(0, 0));
		assertEquals(c.isOutlier(1, 0), false);
		assertTrue(c.isOutlier(5, 0));
		assertTrue(c.isOutlier(0));
		assertEquals(c.isOutlier(1), false);
		assertTrue(c.isOutlier(5));
	}

	/**
	 * Test of getEntry method, of class affymetrix.gcos.cel.FusionCELData.
	 */
	public void testASCII() throws Exception {

		FusionCELData c = new FusionCELData();
		c.setFileName(new File("..\\..\\..\\CPPTest\\data\\test.v3.CEL").getCanonicalPath());
		assertTrue(c.read());
		checkHeader(c, true);
		checkBody(c);
	}

	/**
	 * Test of exists method, of class affymetrix.gcos.cel.FusionCELData.
	 */
	public void testExists() throws Exception {
		FusionCELData c = new FusionCELData();
		c.setFileName(new File("..\\..\\..\\CPPTest\\data\\test.v3.CEL").getCanonicalPath());
		assertTrue(c.exists());
		c.setFileName(new File("..\\..\\..\\CPPTest\\data\\no_file.CEL").getCanonicalPath());
		assertEquals(c.exists(), false);
	}

	/**
	 * Test of readHeader method, of class affymetrix.gcos.cel.FusionCELData.
	 */
	public void testReadHeader() throws Exception {
		FusionCELData c = new FusionCELData();
		c.setFileName(new File("..\\..\\..\\CPPTest\\data\\test.v3.CEL").getCanonicalPath());
		assertTrue(c.readHeader());
		checkHeader(c, false);

		c = new FusionCELData();
		c.setFileName(new File("..\\..\\..\\CPPTest\\data\\test.v4.CEL").getCanonicalPath());
		assertTrue(c.readHeader());
		checkHeader(c, true);

	}

	/**
	 * Test of read method, of class affymetrix.gcos.cel.FusionCELData.
	 */
	public void testXDA() throws Exception {

		FusionCELData c = new FusionCELData();
		c.setFileName(new File("..\\..\\..\\CPPTest\\data\\test.v4.CEL").getCanonicalPath());
		assertTrue(c.read());
		checkHeader(c, true);
		checkBody(c);
	}

}
