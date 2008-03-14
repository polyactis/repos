/*
 * ProbeSetMultiDataCopyNumberDataTest.java
 * JUnit based test
 *
 * Created on June 6, 2007, 12:28 PM
 */

package affymetrix.calvin.data;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * 
 * @author ljevon
 */
public class ProbeSetMultiDataCopyNumberDataTest extends TestCase {

	public ProbeSetMultiDataCopyNumberDataTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(ProbeSetMultiDataCopyNumberDataTest.class);

		return suite;
	}

	public void testName() {
		ProbeSetMultiDataCopyNumberData cnd = new ProbeSetMultiDataCopyNumberData();
		cnd.setName("test");
		assertEquals(cnd.getName(), "test");
	}

	public void testChr() {
		ProbeSetMultiDataCopyNumberData cnd = new ProbeSetMultiDataCopyNumberData();
		cnd.setChr((byte)12);
		assertEquals(cnd.getChr(), 12);
	}

	public void testPosition() {
		ProbeSetMultiDataCopyNumberData cnd = new ProbeSetMultiDataCopyNumberData();
		cnd.setPosition(123);
		assertEquals(cnd.getPosition(), 123);
	}

	public void testChromosomeFromString() {
		assertEquals(ProbeSetMultiDataCopyNumberData.chromosomeFromString("X"), ProbeSetMultiDataCopyNumberData.X_CHR);
		assertEquals(ProbeSetMultiDataCopyNumberData.chromosomeFromString("Y"), ProbeSetMultiDataCopyNumberData.Y_CHR);
		assertEquals(ProbeSetMultiDataCopyNumberData.chromosomeFromString("MT"), ProbeSetMultiDataCopyNumberData.MT_CHR);
		assertEquals(ProbeSetMultiDataCopyNumberData.chromosomeFromString("--"), ProbeSetMultiDataCopyNumberData.NO_CHR);
		assertEquals(ProbeSetMultiDataCopyNumberData.chromosomeFromString("-"), ProbeSetMultiDataCopyNumberData.NO_CHR);
		assertEquals(ProbeSetMultiDataCopyNumberData.chromosomeFromString(""), ProbeSetMultiDataCopyNumberData.NO_CHR);
		assertEquals(ProbeSetMultiDataCopyNumberData.chromosomeFromString("1"), 1);
		assertEquals(ProbeSetMultiDataCopyNumberData.chromosomeFromString("2"), 2);
		assertEquals(ProbeSetMultiDataCopyNumberData.chromosomeFromString("3"), 3);
		assertEquals(ProbeSetMultiDataCopyNumberData.chromosomeFromString("4"), 4);
		assertEquals(ProbeSetMultiDataCopyNumberData.chromosomeFromString("5"), 5);
		assertEquals(ProbeSetMultiDataCopyNumberData.chromosomeFromString("6"), 6);
		assertEquals(ProbeSetMultiDataCopyNumberData.chromosomeFromString("7"), 7);
		assertEquals(ProbeSetMultiDataCopyNumberData.chromosomeFromString("8"), 8);
		assertEquals(ProbeSetMultiDataCopyNumberData.chromosomeFromString("9"), 9);
		assertEquals(ProbeSetMultiDataCopyNumberData.chromosomeFromString("10"), 10);
	}

	public void testChromosomeToString() {
		assertEquals(ProbeSetMultiDataCopyNumberData.chromosomeToString(ProbeSetMultiDataCopyNumberData.X_CHR), "X");
		assertEquals(ProbeSetMultiDataCopyNumberData.chromosomeToString(ProbeSetMultiDataCopyNumberData.Y_CHR), "Y");
		assertEquals(ProbeSetMultiDataCopyNumberData.chromosomeToString(ProbeSetMultiDataCopyNumberData.MT_CHR), "MT");
		assertEquals(ProbeSetMultiDataCopyNumberData.chromosomeToString(ProbeSetMultiDataCopyNumberData.NO_CHR), "-");
		assertEquals(ProbeSetMultiDataCopyNumberData.chromosomeToString((byte)1), "1");
		assertEquals(ProbeSetMultiDataCopyNumberData.chromosomeToString((byte)2), "2");
		assertEquals(ProbeSetMultiDataCopyNumberData.chromosomeToString((byte)3), "3");
		assertEquals(ProbeSetMultiDataCopyNumberData.chromosomeToString((byte)4), "4");
		assertEquals(ProbeSetMultiDataCopyNumberData.chromosomeToString((byte)5), "5");
		assertEquals(ProbeSetMultiDataCopyNumberData.chromosomeToString((byte)6), "6");
		assertEquals(ProbeSetMultiDataCopyNumberData.chromosomeToString((byte)7), "7");
		assertEquals(ProbeSetMultiDataCopyNumberData.chromosomeToString((byte)8), "8");
		assertEquals(ProbeSetMultiDataCopyNumberData.chromosomeToString((byte)9), "9");
		assertEquals(ProbeSetMultiDataCopyNumberData.chromosomeToString((byte)10), "10");
	}

}
