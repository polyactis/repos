/*
 * XYCoordTest.java
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
public class XYCoordTest extends TestCase {

	public XYCoordTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(XYCoordTest.class);

		return suite;
	}

	public void testConstructor() {
		XYCoord xy = new XYCoord();
		assertEquals(xy.getX(), 0);
		assertEquals(xy.getY(), 0);
		xy = new XYCoord((short)1, (short)2);
		assertEquals(xy.getX(), 1);
		assertEquals(xy.getY(), 2);
		XYCoord xy2 = new XYCoord(xy);
		assertEquals(xy2.getX(), 1);
		assertEquals(xy2.getY(), 2);
	}

	public void testProperties() {
		XYCoord xy = new XYCoord();
		xy.setX((short)2);
		assertEquals(xy.getX(), 2);
		xy.setY((short)3);
		assertEquals(xy.getY(), 3);
	}

	public void testLessThan() {
		XYCoord xy1 = new XYCoord((short)1, (short)2);
		XYCoord xy2 = new XYCoord((short)2, (short)4);
		assertTrue(xy1.lessThan(xy2) == true);
		assertTrue(xy2.lessThan(xy1) == false);
	}

}
