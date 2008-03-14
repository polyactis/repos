/*
 * PointU16_tTest.java
 * JUnit based test
 *
 * Created on December 5, 2005, 9:27 AM
 */

package affymetrix.calvin.utils;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * 
 * @author ljevon
 */
public class PointU16_tTest extends TestCase {

	public PointU16_tTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(PointU16_tTest.class);

		return suite;
	}

	/**
	 * Test of getY method, of class affymetrix.calvin.utils.FPoint.
	 */
	public void testProperties() {
		PointU16_t pt = new PointU16_t();
		pt.setX((short)1);
		assertTrue(pt.getX() == 1);
		pt.setY((short)2);
		assertTrue(pt.getY() == 2);
	}

	/**
	 * Test of equals method, of class affymetrix.calvin.utils.FPoint.
	 */
	public void testEquals() {
		PointU16_t pt1 = new PointU16_t();
		pt1.setX((short)2);
		pt1.setY((short)3);
		PointU16_t pt2 = new PointU16_t(pt1);
		assertTrue(pt1.equals(pt2) == true);

	}
}
