/*
 * GridCoordsTest.java
 * JUnit based test
 *
 * Created on December 5, 2005, 9:02 AM
 */

package affymetrix.calvin.utils;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * 
 * @author ljevon
 */
public class GridCoordsTest extends TestCase {

	public GridCoordsTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(GridCoordsTest.class);

		return suite;
	}

	/**
	 * Test of isEmpty method, of class affymetrix.calvin.utils.GridCoords.
	 */
	public void testIsEmpty() {
		GridCoords grid = new GridCoords();
		assertEquals(grid.isEmpty(), true);
		Point pt = new Point();
		pt.setX(1);
		pt.setY(1);
		grid.setLowerLeft(pt);
		grid.setLowerRight(pt);
		grid.setUpperRight(pt);
		grid.setUpperLeft(pt);
		grid.setLowerLeft(pt);
		assertEquals(grid.isEmpty(), true);
	}

	public void testRegion() {
		GridCoords grid = new GridCoords();
		Point pt = new Point();
		pt.setX(1);
		pt.setY(2);
		grid.setUpperLeft(pt);
		pt = new Point();
		pt.setX(3);
		pt.setY(4);
		grid.setUpperRight(pt);
		pt = new Point();
		pt.setX(5);
		pt.setY(6);
		grid.setLowerRight(pt);
		pt = new Point();
		pt.setX(7);
		pt.setY(8);
		grid.setLowerLeft(pt);
		assertEquals(grid.getUpperLeft().getX(), 1);
		assertEquals(grid.getUpperLeft().getY(), 2);
		assertEquals(grid.getUpperRight().getX(), 3);
		assertEquals(grid.getUpperRight().getY(), 4);
		assertEquals(grid.getLowerRight().getX(), 5);
		assertEquals(grid.getLowerRight().getY(), 6);
		assertEquals(grid.getLowerLeft().getX(), 7);
		assertEquals(grid.getLowerLeft().getY(), 8);
		Region r = grid.getRegion();
		assertEquals(r.size(), 4);
		assertEquals(r.get(0).getX(), 1);
		assertEquals(r.get(0).getY(), 2);
		assertEquals(r.get(1).getX(), 3);
		assertEquals(r.get(1).getY(), 4);
		assertEquals(r.get(2).getX(), 5);
		assertEquals(r.get(2).getY(), 6);
		assertEquals(r.get(3).getX(), 7);
		assertEquals(r.get(3).getY(), 8);

	}

}
