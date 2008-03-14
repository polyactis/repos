/*
 * FGridCoordsTest.java
 * JUnit based test
 *
 * Created on November 7, 2005, 9:49 AM
 */

package affymetrix.calvin.utils;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * 
 * @author ljevon
 */
public class FGridCoordsTest extends TestCase {

	public FGridCoordsTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(FGridCoordsTest.class);

		return suite;
	}

	/**
	 * Test of isEmpty method, of class affymetrix.calvin.utils.FGridCoords.
	 */
	public void testIsEmpty() {
		FGridCoords grid = new FGridCoords();
		assertEquals(grid.isEmpty(), true);
		FPoint pt = new FPoint();
		pt.setX(1.0f);
		pt.setY(1.0f);
		grid.setLowerLeft(pt);
		grid.setLowerRight(pt);
		grid.setUpperRight(pt);
		grid.setUpperLeft(pt);
		grid.setLowerLeft(pt);
		assertEquals(grid.isEmpty(), true);
	}

	public void testRegion() {
		FGridCoords grid = new FGridCoords();
		FPoint pt = new FPoint();
		pt.setX(1.0f);
		pt.setY(2.0f);
		grid.setUpperLeft(pt);
		pt = new FPoint();
		pt.setX(3.0f);
		pt.setY(4.0f);
		grid.setUpperRight(pt);
		pt = new FPoint();
		pt.setX(5.0f);
		pt.setY(6.0f);
		grid.setLowerRight(pt);
		pt = new FPoint();
		pt.setX(7.0f);
		pt.setY(8.0f);
		grid.setLowerLeft(pt);

		assertEquals(grid.getUpperLeft().getX(), 1.0f, 0.001);
		assertEquals(grid.getUpperLeft().getY(), 2.0f, 0.001);

		assertEquals(grid.getUpperRight().getX(), 3.0f, 0.001);
		assertEquals(grid.getUpperRight().getY(), 4.0f, 0.001);

		assertEquals(grid.getLowerRight().getX(), 5.0f, 0.001);
		assertEquals(grid.getLowerRight().getY(), 6.0f, 0.001);

		assertEquals(grid.getLowerLeft().getX(), 7.0f, 0.001);
		assertEquals(grid.getLowerLeft().getY(), 8.0f, 0.001);

		FRegion r = grid.getRegion();
		assertEquals(r.size(), 4);
		assertEquals(r.get(0).getX(), 1.0f, 0.001f);
		assertEquals(r.get(0).getY(), 2.0f, 0.001f);
		assertEquals(r.get(1).getX(), 3.0f, 0.001f);
		assertEquals(r.get(1).getY(), 4.0f, 0.001f);
		assertEquals(r.get(2).getX(), 5.0f, 0.001f);
		assertEquals(r.get(2).getY(), 6.0f, 0.001f);
		assertEquals(r.get(3).getX(), 7.0f, 0.001f);
		assertEquals(r.get(3).getY(), 8.0f, 0.001f);

	}

}
