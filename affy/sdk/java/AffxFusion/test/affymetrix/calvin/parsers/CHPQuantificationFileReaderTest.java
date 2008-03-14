/*
 * CHPQuantificationFileReaderTest.java
 * JUnit based test
 *
 * Created on December 27, 2005, 8:21 AM
 */

package affymetrix.calvin.parsers;

import java.io.File;
import java.io.IOException;
import java.util.List;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import affymetrix.calvin.data.CHPQuantificationData;
import affymetrix.calvin.data.ProbeSetQuantificationData;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.calvin.parameter.ParameterNameValue;

/**
 * 
 * @author ljevon
 */
public class CHPQuantificationFileReaderTest extends TestCase {

	public CHPQuantificationFileReaderTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(CHPQuantificationFileReaderTest.class);

		return suite;
	}

	/**
	 * Test of read method, of class affymetrix.calvin.parsers.CHPQuantificationFileReader.
	 */
	public void testRead() throws Exception {
		CHPQuantificationData data = new CHPQuantificationData();
		CHPQuantificationFileReader reader = new CHPQuantificationFileReader();
		String filename = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\CHP_quantification_file")
				.getCanonicalPath();
		reader.setFilename(filename);
		try {
			reader.read(data);
		}
		catch (Throwable t) {
			assertTrue(false);
		}
		assertEquals(data.getFilename(), filename);
		assertEquals(data.getAlgName(), "sig");
		assertEquals(data.getAlgVersion(), "1.0");
		assertEquals(data.getArrayType(), "test3");
		assertEquals(data.getEntryCount(), 2);

		List<ParameterNameValue> params = data.getAlgParams();
		assertEquals(params.size(), 1);
		ParameterNameValue param = params.get(0);
		assertEquals(param.getName(), "an1");
		assertEquals(param.getValueText(), "av1");

		params = data.getSummaryParams();
		assertEquals(params.size(), 1);
		param = params.get(0);
		assertEquals(param.getName(), "sn1");
		assertEquals(param.getValueText(), "sv1");

		ProbeSetQuantificationData e = new ProbeSetQuantificationData();
		e = data.getQuantificationEntry(0);
		assertEquals(e.getName(), "abc");
		assertEquals(e.getId(), -1);
		assertEquals(e.getQuantification(), 10.0f, 0.0001f);
		e = data.getQuantificationEntry(1);
		assertEquals(e.getName(), "xyz");
		assertEquals(e.getId(), -1);
		assertEquals(e.getQuantification(), 20.0f, 0.0001f);

	}

	/**
	 * Test of read method, of class affymetrix.calvin.parsers.CHPQuantificationFileReader.
	 */
	public void testReadId() throws IOException, UnsignedOutOfLimitsException {
		CHPQuantificationData data = new CHPQuantificationData();
		CHPQuantificationFileReader reader = new CHPQuantificationFileReader();
		String filename = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\CHP_quantification_file_id")
				.getCanonicalPath();
		reader.setFilename(filename);
		try {
			reader.read(data);
		}
		catch (Throwable t) {
			assertEquals(true, false);
		}
		assertEquals(data.getFilename(), filename);
		assertEquals(data.getAlgName(), "sig");
		assertEquals(data.getAlgVersion(), "1.0");
		assertEquals(data.getArrayType(), "test3");
		assertEquals(data.getEntryCount(), 2);

		List<ParameterNameValue> params = data.getAlgParams();
		assertEquals(params.size(), 1);
		ParameterNameValue param = params.get(0);
		assertEquals(param.getName(), "an1");
		assertEquals(param.getValueText(), "av1");

		params = data.getSummaryParams();
		assertEquals(params.size(), 1);
		param = params.get(0);
		assertEquals(param.getName(), "sn1");
		assertEquals(param.getValueText(), "sv1");

		ProbeSetQuantificationData e = new ProbeSetQuantificationData();
		e = data.getQuantificationEntry(0);
		assertEquals(e.getName(), "");
		assertEquals(e.getId(), 10);
		assertEquals(e.getQuantification(), 10.0f, 0.0001f);
		e = data.getQuantificationEntry(1);
		assertEquals(e.getName(), "");
		assertEquals(e.getId(), 20);
		assertEquals(e.getQuantification(), 20.0f, 0.0001f);

	}

}
