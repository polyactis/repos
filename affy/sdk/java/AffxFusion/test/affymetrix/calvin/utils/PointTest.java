/*
 * PointTest.java
 * JUnit based test
 *
 * Created on December 5, 2005, 9:25 AM
 */

package affymetrix.calvin.utils;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * 
 * @author ljevon
 */
public class PointTest extends TestCase {

	public PointTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(PointTest.class);

		return suite;
	}

	/**
	 * Test of getY method, of class affymetrix.calvin.utils.FPoint.
	 */
	public void testProperties() {
		Point pt = new Point();
		pt.setX(1);
		assertTrue(pt.getX() == 1);
		pt.setY(2);
		assertTrue(pt.getY() == 2);
	}

	/**
	 * Test of equals method, of class affymetrix.calvin.utils.FPoint.
	 */
	public void testEquals() {
		Point pt1 = new Point();
		pt1.setX(2);
		pt1.setY(3);
		Point pt2 = new Point(pt1);
		assertTrue(pt1.equals(pt2) == true);

	}

}
