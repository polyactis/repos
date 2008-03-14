/*
 * CELFileDataTest.java
 * JUnit based test
 *
 * Created on October 17, 2005, 4:55 PM
 */

package affymetrix.gcos.cel;

import java.io.File;
import java.util.List;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import affymetrix.gcos.TagValuePair;

/**
 * 
 * @author ljevon
 */
public class CELFileDataTest extends TestCase {

	public CELFileDataTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(CELFileDataTest.class);

		return suite;
	}

	/**
	 * Test of getError method, of class affymetrix.gcos.cel.CELFileData.
	 */
	public void testGetError() {
	}

	/**
	 * Test of getFileName method, of class affymetrix.gcos.cel.CELFileData.
	 */
	public void testFileName() {
		CELFileData c = new CELFileData();
		c.setFileName("asdf");
		assertEquals(c.getFileName(), "asdf");
	}

	private void checkHeader(CELFileData c, boolean checkMaskOutliers) {
		assertEquals(c.getHeader().getRows(), 126);
		assertEquals(c.getHeader().getCols(), 126);
		assertEquals(c.getHeader().getCells(), 126 * 126);
		if (checkMaskOutliers == true) {
			assertEquals(c.getHeader().getMasked(), 1);
			assertEquals(c.getHeader().getOutliers(), 3917);
		}

		List<TagValuePair> p = c.getHeader().getParameters();
		String[] ptags = { "Percentile", "CellMargin", "OutlierHigh", "OutlierLow", "AlgVersion", "FixedCellSize",
				"IgnoreOutliersInShiftRows", "FeatureExtraction", "UseSubgrids", "RandomizePixels", "ErrorBasis", "StdMult" };
		String[] pvalues = { "75", "2", "1.500", "1.004", "6.0", "FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "StdvMean",
				"1.000000" };

		assertEquals(p.size(), pvalues.length);
		for (int i = 0; i < p.size(); i++) {
			TagValuePair param = p.get(i);
			assertEquals(param.getTag(), ptags[i]);
			assertEquals(param.getValue(), pvalues[i]);
		}

		assertEquals(c.getHeader().getAlg(), "Percentile");
		assertEquals(c.getHeader().getChipType(), "Test3");
		assertEquals(c.getHeader().getMargin(), 2);
		assertEquals(c.getHeader().getGrid().getUpperLeft().getX(), 239);
		assertEquals(c.getHeader().getGrid().getUpperLeft().getY(), 234);
		assertEquals(c.getHeader().getGrid().getUpperRight().getX(), 4504);
		assertEquals(c.getHeader().getGrid().getUpperRight().getY(), 232);
		assertEquals(c.getHeader().getGrid().getLowerRight().getX(), 4505);
		assertEquals(c.getHeader().getGrid().getLowerRight().getY(), 4497);
		assertEquals(c.getHeader().getGrid().getLowerLeft().getX(), 240);
		assertEquals(c.getHeader().getGrid().getLowerLeft().getY(), 4500);
		assertEquals(
				c.getHeader().getDatHeader(),
				"[81..46222]  Test3:CLS=4733 RWS=4733 XIN=3  YIN=3  VE=17        2.0 02/25/ 0 10:35:25       GridVerify=None  Test3.1sq                  6");
	}

	private void checkBody(CELFileData c) {
		CELFileEntryType entry = new CELFileEntryType();
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
		assertEquals(CELFileData.xyToIndex(0, 0, 126, 126), 0);
		assertEquals(CELFileData.xyToIndex(1, 1, 126, 126), 127);
		assertEquals(CELFileData.xyToIndex(1, 0, 126, 126), 1);

		assertEquals(c.isMasked(0), true);
		assertEquals(c.isMasked(1), false);
		assertEquals(c.isMasked(0, 0), true);
		assertEquals(c.isMasked(1, 0), false);

		assertEquals(c.isOutlier(0, 0), true);
		assertEquals(c.isOutlier(1, 0), false);
		assertEquals(c.isOutlier(5, 0), true);
		assertEquals(c.isOutlier(0), true);
		assertEquals(c.isOutlier(1), false);
		assertEquals(c.isOutlier(5), true);
	}

	/**
	 * Test of getEntry method, of class affymetrix.gcos.cel.CELFileData.
	 */
	public void testASCII() throws Exception {
		CELFileData c = new CELFileData();
		c.setFileName(new File("..\\..\\..\\CPPTest\\data\\test.v3.CEL").getCanonicalPath());
		assertEquals(c.read(), true);
		checkHeader(c, true);
		checkBody(c);
	}

	/**
	 * Test of exists method, of class affymetrix.gcos.cel.CELFileData.
	 */
	public void testExists() throws Exception {
		CELFileData c = new CELFileData();
		c.setFileName(new File("..\\..\\..\\CPPTest\\data\\test.v3.CEL").getCanonicalPath());
		assertEquals(c.exists(), true);
		c.setFileName(new File("..\\..\\..\\CPPTest\\data\\no_file.CEL").getCanonicalPath());
		assertEquals(c.exists(), false);
	}

	/**
	 * Test of readHeader method, of class affymetrix.gcos.cel.CELFileData.
	 */
	public void testReadHeader() throws Exception {
		CELFileData c = new CELFileData();
		c.setFileName(new File("..\\..\\..\\CPPTest\\data\\test.v3.CEL").getCanonicalPath());
		assertEquals(c.readHeader(), true);
		checkHeader(c, false);

		c = new CELFileData();
		c.setFileName(new File("..\\..\\..\\CPPTest\\data\\test.v4.CEL").getCanonicalPath());
		assertEquals(c.readHeader(), true);
		checkHeader(c, true);

	}

	/**
	 * Test of read method, of class affymetrix.gcos.cel.CELFileData.
	 */
	public void testXDA() throws Exception {

		CELFileData c = new CELFileData();
		c.setFileName(new File("..\\..\\..\\CPPTest\\data\\test.v4.CEL").getCanonicalPath());
		assertEquals(c.read(), true);
		checkHeader(c, true);
		checkBody(c);
	}

	/**
	 * Test of isXDACompatibleFile method, of class affymetrix.gcos.cel.CELFileData.
	 */
	public void testIsXDACompatibleFile() throws Exception {
		CELFileData c = new CELFileData();
		c.setFileName(new File("..\\..\\..\\CPPTest\\data\\test.v3.CEL").getCanonicalPath());
		assertEquals(c.isXDACompatibleFile(), false);
		c.setFileName(new File("..\\..\\..\\CPPTest\\data\\test.v4.CEL").getCanonicalPath());
		assertEquals(c.isXDACompatibleFile(), true);
	}

	/**
	 * Test of isVersion3CompatibleFile method, of class affymetrix.gcos.cel.CELFileData.
	 */
	public void testIsVersion3CompatibleFile() throws Exception {
		CELFileData c = new CELFileData();
		c.setFileName(new File("..\\..\\..\\CPPTest\\data\\test.v3.CEL").getCanonicalPath());
		assertEquals(c.isVersion3CompatibleFile(), true);
		c.setFileName(new File("..\\..\\..\\CPPTest\\data\\test.v4.CEL").getCanonicalPath());
		assertEquals(c.isVersion3CompatibleFile(), false);
	}

}
