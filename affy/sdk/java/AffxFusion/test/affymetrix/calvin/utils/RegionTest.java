/*
 * RegionTest.java
 * JUnit based test
 *
 * Created on December 5, 2005, 9:26 AM
 */

package affymetrix.calvin.utils;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * 
 * @author ljevon
 */
public class RegionTest extends TestCase {

	public RegionTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(RegionTest.class);

		return suite;
	}

	/**
	 * Test of clear method, of class affymetrix.calvin.utils.Region.
	 */
	public void testClear() {
		Region r = new Region();
		Point pt = new Point();
		pt.setX(1);
		pt.setY(2);
		r.add(pt);
		assertEquals(r.size(), 1);
		r.clear();
		assertEquals(r.size(), 0);
		r.add(pt);
		pt = null;
		pt = r.get(0);
		assertEquals(pt.getX(), 1.0, 0.001);
		assertEquals(pt.getY(), 2.0, 0.001);
	}

	/**
	 * Test of equals method, of class affymetrix.calvin.utils.Region.
	 */
	public void testEquals() {
		Region r = new Region();
		Point pt = new Point();
		pt.setX(1);
		pt.setY(2);
		r.add(pt);
		Region r2 = new Region();
		r2.add(pt);
		assertTrue(r.equals(r2) == true);

	}

}
