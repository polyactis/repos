/*
 * CELFileHeaderDataTest.java
 * JUnit based test
 *
 * Created on October 17, 2005, 4:55 PM
 */

package affymetrix.gcos.cel;

import java.util.List;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import affymetrix.gcos.CoordinatePoint;
import affymetrix.gcos.GridCoordinates;
import affymetrix.gcos.TagValuePair;

/**
 * 
 * @author ljevon
 */
public class CELFileHeaderDataTest extends TestCase {

	public CELFileHeaderDataTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(CELFileHeaderDataTest.class);

		return suite;
	}

	/**
	 * Test of getMagic method, of class affymetrix.gcos.cel.CELFileHeaderData.
	 */
	public void testMagic() {
		CELFileHeaderData h = new CELFileHeaderData();
		h.setMagic(10);
		assertEquals(h.getMagic(), 10);
	}

	/**
	 * Test of getVersion method, of class affymetrix.gcos.cel.CELFileHeaderData.
	 */
	public void testVersion() {

		CELFileHeaderData h = new CELFileHeaderData();
		h.setVersion(10);
		assertEquals(h.getVersion(), 10);
	}

	/**
	 * Test of getCols method, of class affymetrix.gcos.cel.CELFileHeaderData.
	 */
	public void testCols() {
		CELFileHeaderData h = new CELFileHeaderData();
		h.setCols(10);
		assertEquals(h.getCols(), 10);
	}

	/**
	 * Test of getRows method, of class affymetrix.gcos.cel.CELFileHeaderData.
	 */
	public void testRows() {

		CELFileHeaderData h = new CELFileHeaderData();
		h.setRows(10);
		assertEquals(h.getRows(), 10);
	}

	/**
	 * Test of getCells method, of class affymetrix.gcos.cel.CELFileHeaderData.
	 */
	public void testCells() {
		CELFileHeaderData h = new CELFileHeaderData();
		h.setCells(10);
		assertEquals(h.getCells(), 10);
	}

	/**
	 * Test of getHeader method, of class affymetrix.gcos.cel.CELFileHeaderData.
	 */
	public void testHeader() {

		CELFileHeaderData h = new CELFileHeaderData();
		h.setHeader("abc");
		assertEquals(h.getHeader(), "abc");
	}

	/**
	 * Test of getAlg method, of class affymetrix.gcos.cel.CELFileHeaderData.
	 */
	public void testAlg() {
		CELFileHeaderData h = new CELFileHeaderData();
		h.setAlg("alg");
		assertEquals(h.getAlg(), "alg");
	}

	/**
	 * Test of getChipType method, of class affymetrix.gcos.cel.CELFileHeaderData.
	 */
	public void testChipType() {
		CELFileHeaderData h = new CELFileHeaderData();
		h.setChipType("test3");
		assertEquals(h.getChipType(), "test3");
	}

	/**
	 * Test of getDatHeader method, of class affymetrix.gcos.cel.CELFileHeaderData.
	 */
	public void testDatHeader() {
		CELFileHeaderData h = new CELFileHeaderData();
		h.setDatHeader("dat");
		assertEquals(h.getDatHeader(), "dat");
	}

	/**
	 * Test of getMargin method, of class affymetrix.gcos.cel.CELFileHeaderData.
	 */
	public void testMargin() {
		CELFileHeaderData h = new CELFileHeaderData();
		h.setMargin(10);
		assertEquals(h.getMargin(), 10);
	}

	/**
	 * Test of getOutliers method, of class affymetrix.gcos.cel.CELFileHeaderData.
	 */
	public void testOutliers() {
		CELFileHeaderData h = new CELFileHeaderData();
		h.setOutliers(10);
		assertEquals(h.getOutliers(), 10);
	}

	/**
	 * Test of getMasked method, of class affymetrix.gcos.cel.CELFileHeaderData.
	 */
	public void testMasked() {
		CELFileHeaderData h = new CELFileHeaderData();
		h.setMasked(10);
		assertEquals(h.getMasked(), 10);
	}

	/**
	 * Test of getGrid method, of class affymetrix.gcos.cel.CELFileHeaderData.
	 */
	public void testGrid() {
		CELFileHeaderData h = new CELFileHeaderData();
		assertEquals(h.getGrid(), null);
		CoordinatePoint pt = new CoordinatePoint();
		pt.setX(1);
		pt.setY(2);
		GridCoordinates grid = new GridCoordinates();
		grid.setUpperRight(new CoordinatePoint(pt));
		pt.setX(2);
		pt.setY(3);
		grid.setUpperLeft(new CoordinatePoint(pt));
		pt.setX(3);
		pt.setY(4);
		grid.setLowerRight(new CoordinatePoint(pt));
		pt.setX(4);
		pt.setY(5);
		grid.setLowerLeft(new CoordinatePoint(pt));
		h.setGrid(grid);
		grid = null;

		grid = h.getGrid();
		assertEquals(grid.getUpperRight().getX(), 1);
		assertEquals(grid.getUpperRight().getY(), 2);
		assertEquals(grid.getUpperLeft().getX(), 2);
		assertEquals(grid.getUpperLeft().getY(), 3);
		assertEquals(grid.getLowerRight().getX(), 3);
		assertEquals(grid.getLowerRight().getY(), 4);
		assertEquals(grid.getLowerLeft().getX(), 4);
		assertEquals(grid.getLowerLeft().getY(), 5);
	}

	/**
	 * Test of getParameters method, of class affymetrix.gcos.cel.CELFileHeaderData.
	 */
	public void testParameters() {
		CELFileHeaderData h = new CELFileHeaderData();
		h.setParameters("A=1 B=2");
		List<TagValuePair> p = h.getParameters();
		TagValuePair param = p.get(0);
		assertEquals(param.getTag(), "A");
		assertEquals(param.getValue(), "1");
		param = p.get(1);
		assertEquals(param.getTag(), "B");
		assertEquals(param.getValue(), "2");
		h.setParameters("A:1;B:3");
		p = h.getParameters();
		param = p.get(0);
		assertEquals(param.getTag(), "A");
		assertEquals(param.getValue(), "1");
		param = p.get(1);
		assertEquals(param.getTag(), "B");
		assertEquals(param.getValue(), "3");
	}

	/**
	 * Test of parseChipType method, of class affymetrix.gcos.cel.CELFileHeaderData.
	 */
	public void testParseChipType() {
		CELFileHeaderData h = new CELFileHeaderData();
		h.setHeader("blah blah blah test3.1sq blah blah blah");
		h.parseChipType();
		assertEquals(h.getChipType(), "test3");
	}

	/**
	 * Test of parseMargin method, of class affymetrix.gcos.cel.CELFileHeaderData.
	 */
	public void testParseMargin() {
		CELFileHeaderData h = new CELFileHeaderData();
		h.setHeader("blah blah blah CellMargin=3 blah blah blah");
		h.parseMargin();
		assertEquals(h.getMargin(), 3);
	}

	/**
	 * Test of parseGrid method, of class affymetrix.gcos.cel.CELFileHeaderData.
	 */
	public void testParseGrid() {
		CELFileHeaderData h = new CELFileHeaderData();
		h.setHeader("blah blah blah GridCornerUL=1 2 GridCornerUR=2 3 GridCornerLR=3 4 GridCornerLL=4 5Axis blah");
		h.parseGrid();
		GridCoordinates grid = h.getGrid();
		assertEquals(grid.getUpperRight().getX(), 2);
		assertEquals(grid.getUpperRight().getY(), 3);
		assertEquals(grid.getUpperLeft().getX(), 1);
		assertEquals(grid.getUpperLeft().getY(), 2);
		assertEquals(grid.getLowerRight().getX(), 3);
		assertEquals(grid.getLowerRight().getY(), 4);
		assertEquals(grid.getLowerLeft().getX(), 4);
		assertEquals(grid.getLowerLeft().getY(), 5);

	}

	/**
	 * Test of parseDatHeader method, of class affymetrix.gcos.cel.CELFileHeaderData.
	 */
	public void testParseDatHeader() {
		CELFileHeaderData h = new CELFileHeaderData();
		h.setHeader("blah blah DatHeader=123312 123321 Algorithm=");
		h.parseDatHeader();
		assertEquals(h.getDatHeader(), "123312 123321");
	}

	/**
	 * Test of clear method, of class affymetrix.gcos.cel.CELFileHeaderData.
	 */
	public void testClear() {
		CELFileHeaderData h = new CELFileHeaderData();
		h.clear();
		assertEquals(h.getAlg(), "");
		assertEquals(h.getCells(), 0);
		assertEquals(h.getRows(), 0);
		assertEquals(h.getCols(), 0);
		assertEquals(h.getGrid(), null);
		assertEquals(h.getChipType(), "");
		assertEquals(h.getDatHeader(), "");
		assertEquals(h.getMargin(), 0);
		assertEquals(h.getMasked(), 0);
		assertEquals(h.getOutliers(), 0);
		assertEquals(h.getParameters().size(), 0);

	}

}
