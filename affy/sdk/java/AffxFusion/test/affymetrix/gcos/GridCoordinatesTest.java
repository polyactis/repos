/*
 * GridCoordinatesTest.java
 * JUnit based test
 *
 * Created on October 16, 2005, 7:55 PM
 */

package affymetrix.gcos;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * 
 * @author ljevon
 */
public class GridCoordinatesTest extends TestCase {

	public GridCoordinatesTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(GridCoordinatesTest.class);

		return suite;
	}

	/**
	 * Test of getUpperLeft method, of class affymetrix.gcos.GridCoordinates.
	 */
	public void testUpperLeft() {
		GridCoordinates grid = new GridCoordinates();
		CoordinatePoint c = new CoordinatePoint();
		c.setX(10);
		c.setY(20);
		grid.setUpperLeft(c);
		c.setX(0);
		c.setY(0);
		CoordinatePoint c1 = grid.getUpperLeft();
		assertEquals(c1.getX(), 10);
		assertEquals(c1.getY(), 20);
	}

	/**
	 * Test of getUppeRright method, of class affymetrix.gcos.GridCoordinates.
	 */
	public void testUpperRight() {
		GridCoordinates grid = new GridCoordinates();
		CoordinatePoint c = new CoordinatePoint();
		c.setX(10);
		c.setY(20);
		grid.setUpperRight(c);
		c.setX(0);
		c.setY(0);
		CoordinatePoint c1 = grid.getUpperRight();
		assertEquals(c1.getX(), 10);
		assertEquals(c1.getY(), 20);
	}

	/**
	 * Test of getLowerRight method, of class affymetrix.gcos.GridCoordinates.
	 */
	public void testLowerRight() {
		GridCoordinates grid = new GridCoordinates();
		CoordinatePoint c = new CoordinatePoint();
		c.setX(10);
		c.setY(20);
		grid.setLowerRight(c);
		c.setX(0);
		c.setY(0);
		CoordinatePoint c1 = grid.getLowerRight();
		assertEquals(c1.getX(), 10);
		assertEquals(c1.getY(), 20);
	}

	/**
	 * Test of getLowerLeft method, of class affymetrix.gcos.GridCoordinates.
	 */
	public void testLowerLeft() {
		GridCoordinates grid = new GridCoordinates();
		CoordinatePoint c = new CoordinatePoint();
		c.setX(10);
		c.setY(20);
		grid.setLowerLeft(c);
		c.setX(0);
		c.setY(0);
		CoordinatePoint c1 = grid.getLowerLeft();
		assertEquals(c1.getX(), 10);
		assertEquals(c1.getY(), 20);
	}

}
