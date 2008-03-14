/*
 * CoordinatePointTest.java
 * JUnit based test
 *
 * Created on October 16, 2005, 7:47 PM
 */

package affymetrix.gcos;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * 
 * @author ljevon
 */
public class CoordinatePointTest extends TestCase {

	public CoordinatePointTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(CoordinatePointTest.class);

		return suite;
	}

	/**
	 * Test of getX method, of class affymetrix.gcos.CoordinatePoint.
	 */
	public void testX() {
		CoordinatePoint c = new CoordinatePoint();
		c.setX(10);
		assertEquals(c.getX(), 10);
		assertEquals(c.getY(), 0);
	}

	/**
	 * Test of getY method, of class affymetrix.gcos.CoordinatePoint.
	 */
	public void testY() {
		CoordinatePoint c = new CoordinatePoint();
		c.setY(10);
		assertEquals(c.getY(), 10);
		assertEquals(c.getX(), 0);
	}
}
