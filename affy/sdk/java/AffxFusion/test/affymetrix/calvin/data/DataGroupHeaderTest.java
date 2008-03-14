/*
 * DataGroupHeaderTest.java
 * JUnit based test
 *
 * Created on October 23, 2005, 7:05 PM
 */

package affymetrix.calvin.data;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import affymetrix.calvin.parameter.ParameterNameValue;

/**
 * 
 * @author ljevon
 */
public class DataGroupHeaderTest extends TestCase {

	public DataGroupHeaderTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(DataGroupHeaderTest.class);

		return suite;
	}

	public void testNameTest() {
		DataGroupHeader header = new DataGroupHeader("dataGroup");
		String p1 = "some_name";
		header.setName(p1);
		String p2 = header.getName();
		assertEquals(p1, p2);
	}

	public void testFindDataSetHeaderTest() {
		// Create DataSetHeaders
		DataGroupHeader header = new DataGroupHeader("dataGroup");
		DataSetHeader dph1 = new DataSetHeader();
		dph1.setName("pixel intensity");
		ParameterNameValue param = new ParameterNameValue();
		param.setName("Scanner");
		param.setValueText("M10");
		dph1.addNameValParam(param);
		dph1.addUShortColumn("");
		dph1.setRowCnt(1);

		DataSetHeader dph2 = new DataSetHeader();
		dph2.setName("grid coordinates");
		param.setName("Corner Pattern");
		param.setValueText("Checkerboard");
		dph2.addNameValParam(param);
		dph2.addUShortColumn("");
		dph2.addUShortColumn("");
		dph2.addUShortColumn("");
		dph2.addUShortColumn("");

		dph2.addUShortColumn("");
		dph2.addUShortColumn("");
		dph2.addUShortColumn("");
		dph2.addUShortColumn("");

		dph2.addUShortColumn("");
		dph2.addUShortColumn("");
		dph2.addUShortColumn("");
		dph2.addUShortColumn("");

		dph2.addUShortColumn("");
		dph2.addUShortColumn("");
		dph2.addUShortColumn("");
		dph2.addUShortColumn("");
		dph2.setRowCnt(1);

		header.addDataSetHdr(dph1);
		header.addDataSetHdr(dph2);

		DataSetHeader dph = header.findDataSetHeader("none");
		assertEquals(dph, null);
		dph = header.findDataSetHeader(dph1.getName());
		assertTrue(null != dph);
		assertEquals(dph.getName(), dph1.getName());
		dph = header.findDataSetHeader(dph2.getName());
		assertTrue(null != dph);
		assertEquals(dph.getName(), dph2.getName());
	}

}
