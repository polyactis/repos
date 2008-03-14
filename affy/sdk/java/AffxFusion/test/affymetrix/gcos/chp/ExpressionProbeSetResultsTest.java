/*
 * ExpressionProbeSetResultsTest.java
 * JUnit based test
 *
 * Created on October 13, 2005, 8:03 AM
 */

package affymetrix.gcos.chp;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import affymetrix.portability.UByte;
import affymetrix.portability.UShort;

/**
 * 
 * @author ljevon
 */
public class ExpressionProbeSetResultsTest extends TestCase {

	public ExpressionProbeSetResultsTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(ExpressionProbeSetResultsTest.class);

		return suite;
	}

	/**
	 * Test of getDetectionPValue method, of class affymetrix.gcos.chp.ExpressionProbeSetResults.
	 */
	public void testDetectionPValue() {
		ExpressionProbeSetResults p = new ExpressionProbeSetResults();
		p.setDetectionPValue(0.5f);
		assertEquals(p.getDetectionPValue(), 0.5f, 0.0000001f);
	}

	/**
	 * Test of getSignal method, of class affymetrix.gcos.chp.ExpressionProbeSetResults.
	 */
	public void testSignal() {
		ExpressionProbeSetResults p = new ExpressionProbeSetResults();
		p.setSignal(0.5f);
		assertEquals(p.getSignal(), 0.5f, 0.0000001f);
	}

	/**
	 * Test of getNumPairs method, of class affymetrix.gcos.chp.ExpressionProbeSetResults.
	 */
	public void testNumPairs() throws Exception {
		ExpressionProbeSetResults p = new ExpressionProbeSetResults();
		p.setNumPairs(new UShort(10));
		assertEquals(p.getNumPairs().toInt(), 10);
	}

	/**
	 * Test of setNumPairs method, of class affymetrix.gcos.chp.ExpressionProbeSetResults.
	 */
	public void testNumUsedPairs() throws Exception {
		ExpressionProbeSetResults p = new ExpressionProbeSetResults();
		p.setNumUsedPairs(new UShort(10));
		assertEquals(p.getNumUsedPairs().toInt(), 10);
	}

	/**
	 * Test of getDetection method, of class affymetrix.gcos.chp.ExpressionProbeSetResults.
	 */
	public void testDetection() throws Exception {
		ExpressionProbeSetResults p = new ExpressionProbeSetResults();
		p.setDetection(new UByte(ExpressionProbeSetResults.ABS_ABSENT_CALL));
		assertEquals(p.getDetection().toShort(), ExpressionProbeSetResults.ABS_ABSENT_CALL);
	}

	/**
	 * Test of getHasCompResults method, of class affymetrix.gcos.chp.ExpressionProbeSetResults.
	 */
	public void testHasCompResults() {
		ExpressionProbeSetResults p = new ExpressionProbeSetResults();
		assertEquals(p.getHasCompResults(), false);
		p.setHasCompResults(true);
		assertEquals(p.getHasCompResults(), true);
	}

	/**
	 * Test of getChangePValue method, of class affymetrix.gcos.chp.ExpressionProbeSetResults.
	 */
	public void testChangePValue() {
		ExpressionProbeSetResults p = new ExpressionProbeSetResults();
		p.setChangePValue(0.5f);
		assertEquals(p.getChangePValue(), 0.5f, 0.0000001f);
	}

	/**
	 * Test of getSignalLogRatio method, of class affymetrix.gcos.chp.ExpressionProbeSetResults.
	 */
	public void testSignalLogRatio() {
		ExpressionProbeSetResults p = new ExpressionProbeSetResults();
		p.setSignalLogRatio(0.5f);
		assertEquals(p.getSignalLogRatio(), 0.5f, 0.0000001f);
	}

	/**
	 * Test of getSignalLogRatioLow method, of class affymetrix.gcos.chp.ExpressionProbeSetResults.
	 */
	public void testSignalLogRatioLow() {
		ExpressionProbeSetResults p = new ExpressionProbeSetResults();
		p.setSignalLogRatioLow(0.5f);
		assertEquals(p.getSignalLogRatioLow(), 0.5f, 0.0000001f);
	}

	/**
	 * Test of getSignalLogRatioHigh method, of class affymetrix.gcos.chp.ExpressionProbeSetResults.
	 */
	public void testSignalLogRatioHigh() {
		ExpressionProbeSetResults p = new ExpressionProbeSetResults();
		p.setSignalLogRatioHigh(0.5f);
		assertEquals(p.getSignalLogRatioHigh(), 0.5f, 0.0000001f);
	}

	/**
	 * Test of getNumCommonPairs method, of class affymetrix.gcos.chp.ExpressionProbeSetResults.
	 */
	public void testNumCommonPairs() throws Exception {
		ExpressionProbeSetResults p = new ExpressionProbeSetResults();
		p.setNumCommonPairs(new UShort(5));
		assertEquals(p.getNumCommonPairs().toInt(), 5);
	}

	/**
	 * Test of getChange method, of class affymetrix.gcos.chp.ExpressionProbeSetResults.
	 */
	public void testChange() throws Exception {
		ExpressionProbeSetResults p = new ExpressionProbeSetResults();
		p.setChange(new UByte(ExpressionProbeSetResults.COMP_DECREASE_CALL));
		assertEquals(p.getChange().toShort(), ExpressionProbeSetResults.COMP_DECREASE_CALL);
	}

	/**
	 * Test of getDetectionString method, of class affymetrix.gcos.chp.ExpressionProbeSetResults.
	 */
	public void testGetDetectionString() throws Exception {
		ExpressionProbeSetResults p = new ExpressionProbeSetResults();
		p.setDetection(new UByte(ExpressionProbeSetResults.ABS_ABSENT_CALL));
		assertEquals(p.getDetectionString(), "A");
		p.setDetection(new UByte(ExpressionProbeSetResults.ABS_PRESENT_CALL));
		assertEquals(p.getDetectionString(), "P");
		p.setDetection(new UByte(ExpressionProbeSetResults.ABS_MARGINAL_CALL));
		assertEquals(p.getDetectionString(), "M");
		p.setDetection(new UByte(ExpressionProbeSetResults.ABS_NO_CALL));
		assertEquals(p.getDetectionString(), "No Call");
	}

	/**
	 * Test of getChangeString method, of class affymetrix.gcos.chp.ExpressionProbeSetResults.
	 */
	public void testGetChangeString() throws Exception {
		ExpressionProbeSetResults p = new ExpressionProbeSetResults();
		p.setChange(new UByte(ExpressionProbeSetResults.COMP_DECREASE_CALL));
		assertEquals(p.getChangeString(), "D");
		p.setChange(new UByte(ExpressionProbeSetResults.COMP_INCREASE_CALL));
		assertEquals(p.getChangeString(), "I");
		p.setChange(new UByte(ExpressionProbeSetResults.COMP_MOD_DECREASE_CALL));
		assertEquals(p.getChangeString(), "MD");
		p.setChange(new UByte(ExpressionProbeSetResults.COMP_MOD_INCREASE_CALL));
		assertEquals(p.getChangeString(), "MI");
		p.setChange(new UByte(ExpressionProbeSetResults.COMP_NO_CHANGE_CALL));
		assertEquals(p.getChangeString(), "NC");
		p.setChange(new UByte(ExpressionProbeSetResults.COMP_NO_CALL));
		assertEquals(p.getChangeString(), "No Call");
	}

}
