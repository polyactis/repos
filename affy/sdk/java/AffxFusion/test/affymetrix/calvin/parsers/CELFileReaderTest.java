/*
 * CELFileReaderTest.java
 * JUnit based test
 *
 * Created on October 31, 2005, 8:41 AM
 */

package affymetrix.calvin.parsers;

import java.io.File;
import java.util.List;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import affymetrix.calvin.data.CELData;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.calvin.parameter.ParameterNameValue;

/**
 * 
 * @author ljevon
 */
public class CELFileReaderTest extends TestCase {

	public CELFileReaderTest(String testName) throws Exception {
		super(testName);
		SMALL_CEL_FILE = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\small_cel_file").getCanonicalPath();

		SMALL_CEL_FILE_NO_OUTLIER_NO_MASK = new File(
				"..\\..\\..\\..\\calvin_files\\parsers\\data\\small_cel_file_no_outlier_no_mask").getCanonicalPath();

		SMALL_CEL_FILE_NO_STDEV = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\small_cel_file_no_stdev")
				.getCanonicalPath();

		LARGE_CEL_FILE = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\large_cel_file").getCanonicalPath();

		SMALL_UINT16_CEL_FILE = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\small_uint16_cel_file")
				.getCanonicalPath();
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(CELFileReaderTest.class);

		return suite;
	}

	/**
	 * Test of getFilename method, of class affymetrix.calvin.parsers.CELFileReader.
	 */
	public void testFilename() {
		CELFileReader reader = new CELFileReader();
		reader.setFilename("file");
		assertEquals(reader.getFilename(), "file");
	}

	/**
	 * Test of read method, of class affymetrix.calvin.parsers.CELFileReader.
	 */
	public void testRead() throws Exception {
		CELData data = new CELData();
		CELFileReader reader = new CELFileReader();
		String fn = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\small_cel_file").getCanonicalPath();
		reader.setFilename(fn);
		try {
			reader.read(data);
		}
		catch (Throwable t) {
			assertTrue(false);
		}

		// Check values
		int cells = data.getNumCells();
		for (int cell = 0; cell < cells; ++cell) {
			assertEquals(data.getIntensity(cell), 100.0f * cell, 0.000001f);
			assertEquals(data.getStdv(cell), .5 * cell, 0.000001f);
			assertEquals(data.getPixels(cell), 25);
		}

	}

	private String SMALL_CEL_FILE;

	private String SMALL_CEL_FILE_NO_OUTLIER_NO_MASK;

	private String SMALL_CEL_FILE_NO_STDEV;

	private String LARGE_CEL_FILE;

	private String SMALL_UINT16_CEL_FILE;

	private static final int[] intenTestValues = { 55683, 4568, 2368, 100 };

	private static final float[] stdevTestValues = { 2.345f, 56.23f, 1.53f, 3.875f };

	public void testReadSmallCelFileNoOutlierNoMask() throws UnsignedOutOfLimitsException {
		CELData data = new CELData();
		CELFileReader reader = new CELFileReader();
		reader.setFilename(SMALL_CEL_FILE_NO_OUTLIER_NO_MASK);
		try {
			reader.read(data);
		}
		catch (Throwable t) {
			assertTrue(false);
		}

		// Read some data
		assertEquals(data.getFileHeader().getFilename(), SMALL_CEL_FILE_NO_OUTLIER_NO_MASK);
		assertEquals(data.getVersion().toShort(), 1);
		assertEquals(data.getArrayType(), "Hg-small");
		assertEquals(data.getAlgorithmName(), "Feature Extraction");
		assertEquals(data.getRows(), 2);
		assertEquals(data.getCols(), 5);
		assertEquals(data.getNumCells(), 10);
		assertTrue(data.hasStdev());
		assertTrue(data.hasNumPixels());

		List<ParameterNameValue> params = data.getAlgorithmParameters();
		assertEquals(params.size(), 2);
		ParameterNameValue param = params.get(0);
		assertEquals(param.getName(), "percentile");
		assertEquals(param.getValueFloat(), 0.75f, 0.000001f);
		param = params.get(1);
		assertEquals(param.getName(), "outlierlow");
		assertEquals(param.getValueFloat(), 1.004f, 0.000001f);

		int cells = data.getNumCells();

		for (int cell = 0; cell < cells; ++cell) {
			assertEquals(data.getIntensity(cell), 100.0f * cell, 0.00001f);
			assertEquals(data.getStdv(cell), .5 * cell, 0.00001f);
			assertEquals(data.getPixels(cell), 25);

			// Check the values
			assertFalse(data.isOutlier(cell));
			assertFalse(data.isMasked(cell));
		}
	}

	public void testReadSmallCelFile() throws Exception {
		CELData data = new CELData();
		CELFileReader reader = new CELFileReader();
		reader.setFilename(SMALL_CEL_FILE);
		try {
			reader.read(data);
		}
		catch (Throwable t) {
			assertTrue(false);
		}

		assertEquals(data.getFileHeader().getFilename(), SMALL_CEL_FILE);
		checkSmallCelFileHeader(data);

		int cells = data.getNumCells();

		// Check values
		for (int cell = 0; cell < cells; ++cell) {
			assertEquals(data.getIntensity(cell), 100 * cell, 0.00001f);
			assertEquals(data.getStdv(cell), .5 * cell, 0.00001f);
			assertEquals(data.getPixels(cell), 25);
			CheckOutlier(cell, data.isOutlier(cell));
			CheckMasked(cell, data.isMasked(cell));
		}

	}

	void checkSmallCelFileHeader(CELData data) throws Exception {
		// Read some data
		assertEquals(data.getVersion().toShort(), 1);
		assertEquals(data.getArrayType(), "Hg-small");
		assertEquals(data.getAlgorithmName(), "Feature Extraction");
		assertEquals(data.getRows(), 5);
		assertEquals(data.getCols(), 5);
		assertEquals(data.getNumCells(), 25);
		assertTrue(data.hasStdev());
		assertTrue(data.hasNumPixels());

		List<ParameterNameValue> params = data.getAlgorithmParameters();
		assertEquals(params.size(), 2);
		ParameterNameValue param = (ParameterNameValue)params.get(0);
		assertEquals(param.getName(), "percentile");
		assertEquals(param.getValueFloat(), 0.75f, 0.000001f);
		param = (ParameterNameValue)params.get(1);
		assertEquals(param.getName(), "outlierlow");
		assertEquals(param.getValueFloat(), 1.004f, 0.000001f);
	}

	public void testReadSmallCelFileNoStdev() throws UnsignedOutOfLimitsException {
		CELData data = new CELData();
		CELFileReader reader = new CELFileReader();
		reader.setFilename(SMALL_CEL_FILE_NO_STDEV);
		try {
			reader.read(data);
		}
		catch (Throwable t) {
			assertTrue(false);
		}

		// Read some data
		assertEquals(data.getFileHeader().getFilename(), SMALL_CEL_FILE_NO_STDEV);
		assertEquals(data.getVersion().toShort(), 1);
		assertEquals(data.getArrayType(), "Hg-small");
		assertEquals(data.getAlgorithmName(), "Feature Extraction");
		assertEquals(data.getRows(), 5);
		assertEquals(data.getCols(), 5);
		assertEquals(data.getNumCells(), 25);
		assertFalse(data.hasStdev());
		assertTrue(data.hasNumPixels());

		List<ParameterNameValue> params = data.getAlgorithmParameters();
		assertEquals(params.size(), 2);
		ParameterNameValue param = (ParameterNameValue)params.get(0);
		assertEquals(param.getName(), "percentile");
		assertEquals(param.getValueFloat(), 0.75f, 0.000001f);
		param = (ParameterNameValue)params.get(1);
		assertEquals(param.getName(), "outlierlow");
		assertEquals(param.getValueFloat(), 1.004f, 0.000001f);

		// float intensity;
		// short numPixels;
		// boolean outlier, masked;
		// float stdev;

		int cells = data.getNumCells();

		for (int cell = 0; cell < cells; ++cell) {
			assertEquals(data.getIntensity(cell), 100 * cell, 0.00001f);
			assertEquals(data.getStdv(cell), 0, 0.000001f);
			assertEquals(data.getPixels(cell), 25);

			CheckOutlier(cell, data.isOutlier(cell));
			CheckMasked(cell, data.isMasked(cell));
		}

		// Check values
		for (int cell = 0; cell < cells; ++cell) {
			assertEquals(data.getIntensity(cell), 100 * cell, 0.000001f);
			assertEquals(data.getPixels(cell), 25);

			CheckOutlier(cell, data.isOutlier(cell));
			CheckMasked(cell, data.isMasked(cell));
		}

		assertEquals(data.getNumOutliers(), 2);
		assertEquals(data.getNumMasked(), 3);
	}

	/*
	 * public void testReadLargeCelFileCheckHeaderTest() { CELData data = new CELData(); CELFileReader reader = new
	 * CELFileReader(); reader.setFilename(LARGE_CEL_FILE); try { reader.read(data); } catch (Throwable t) {
	 * assertTrue(false); }
	 * 
	 * int cols = 2560; int rows = 2560; int expectedCells = rows*cols; // Read some data
	 * assertEquals(data.getFileHeader().getFilename(), LARGE_CEL_FILE); assertEquals(data.getRows(), rows);
	 * assertEquals(data.getCols(), cols); assertEquals(data.getNumCells(), expectedCells); assertEquals(data.hasStdev(),
	 * true); assertEquals(data.hasNumPixels(), true); }
	 */

	public void testReadLargeCelFileSingleDataAccess() throws UnsignedOutOfLimitsException {
		CELData data = new CELData();
		CELFileReader reader = new CELFileReader();
		reader.setFilename(LARGE_CEL_FILE);
		try {
			reader.read(data);
		}
		catch (Throwable t) {
			assertTrue(false);
		}

		int cols = 2560;
		int rows = 2560;
		int cell = 0;

		// This needs to be faster. Spot check
		int row;
		int col;
		int testIdx;
		for (row = 0, testIdx = 0; row < rows; row += 1) {
			for (col = 0, testIdx = 0; col < cols; ++col, ++testIdx) {
				// Compute the expected values
				if (testIdx >= 4)
					testIdx = 0;

				assertEquals(data.getIntensity(cell), intenTestValues[testIdx], 0.00001f);
				assertEquals(data.getStdv(cell), stdevTestValues[testIdx], 0.00001f);
				assertEquals(data.getPixels(cell), 25);
				++cell;
			}
		}
	}

	void CheckOutlier(int cell, boolean outlier) {
		if (cell == 0 || cell == 11)
			assertTrue(outlier);
		else
			assertFalse(outlier);
	}

	void CheckMasked(int cell, boolean masked) {
		if (cell == 1 || cell == 7 || cell == 13)
			assertTrue(masked);
		else
			assertFalse(masked);
	}

	// Test added to reproduce bug reported by beta testers.
	public void testReadSmallCelFileGetOutlierCoordsTest() {
		CELData data = new CELData();
		CELFileReader reader = new CELFileReader();
		reader.setFilename(SMALL_CEL_FILE);
		try {
			reader.read(data);
		}
		catch (Throwable t) {
			assertTrue(false);
		}

		assertEquals(data.getNumOutliers(), 2);
	}

	// Test added to reproduce bug reported by beta testers.
	public void testReadSmallCelFileGetMaskedCoords() {
		CELData data = new CELData();
		CELFileReader reader = new CELFileReader();
		reader.setFilename(SMALL_CEL_FILE);
		try {
			reader.read(data);
		}
		catch (Throwable t) {
			assertTrue(false);
		}
		assertEquals(data.getNumMasked(), 3);
	}

	public void testReadSmallUInt16CelFile() throws Exception {
		CELData data = new CELData();
		CELFileReader reader = new CELFileReader();
		reader.setFilename(SMALL_UINT16_CEL_FILE);
		try {
			reader.read(data);
		}
		catch (Throwable t) {
			assertTrue(false);
		}

		// Read some data
		assertEquals(data.getFileHeader().getFilename(), SMALL_UINT16_CEL_FILE);

		checkSmallCelFileHeader(data);

		// float intensity;
		// short numPixels;
		// boolean outlier, masked;
		// float stdev;

		int cells = data.getNumCells();

		// Check values
		for (int cell = 0; cell < cells; ++cell) {
			assertEquals(data.getIntensity(cell), 100 * cell, 0.000001f);
			assertEquals(data.getStdv(cell), .5 * cell, 0.00001f);
			assertEquals(data.getPixels(cell), 25);
			CheckOutlier(cell, data.isOutlier(cell));
			CheckMasked(cell, data.isMasked(cell));
		}

		assertEquals(data.getNumOutliers(), 2);
		assertEquals(data.getNumMasked(), 3);

		// Change the order of operations in the test
		for (int cell = 0; cell < cells; ++cell) {
			assertEquals(data.getIntensity(cell), 100 * cell, 0.00001f);
			assertEquals(data.getStdv(cell), .5 * cell, 0.000001f);
			assertEquals(data.getPixels(cell), 25);

			CheckOutlier(cell, data.isOutlier(cell));
			CheckMasked(cell, data.isMasked(cell));
		}
	}

}
