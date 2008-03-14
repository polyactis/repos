/*
 * CHPFileHeaderTest.java
 * JUnit based test
 *
 * Created on October 13, 2005, 8:19 AM
 */

package affymetrix.gcos.chp;

import java.util.ArrayList;
import java.util.List;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import affymetrix.gcos.TagValuePair;

/**
 * 
 * @author ljevon
 */
public class CHPFileHeaderTest extends TestCase {

	public CHPFileHeaderTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(CHPFileHeaderTest.class);
		return suite;
	}

	/**
	 * Test of getMagic method, of class affymetrix.gcos.chp.CHPFileHeader.
	 */
	public void testMagic() {
		CHPFileHeader h = new CHPFileHeader();
		h.setMagic(10);
		assertEquals(h.getMagic(), 10);
	}

	/**
	 * Test of getVersion method, of class affymetrix.gcos.chp.CHPFileHeader.
	 */
	public void testVersion() {
		CHPFileHeader h = new CHPFileHeader();
		h.setVersion(10);
		assertEquals(h.getVersion(), 10);
	}

	/**
	 * Test of getCols method, of class affymetrix.gcos.chp.CHPFileHeader.
	 */
	public void testCols() {
		CHPFileHeader h = new CHPFileHeader();
		h.setCols(10);
		assertEquals(h.getCols(), 10);
	}

	/**
	 * Test of getRows method, of class affymetrix.gcos.chp.CHPFileHeader.
	 */
	public void testRows() {
		CHPFileHeader h = new CHPFileHeader();
		h.setRows(10);
		assertEquals(h.getRows(), 10);
	}

	/**
	 * Test of getNumProbeSets method, of class affymetrix.gcos.chp.CHPFileHeader.
	 */
	public void testNumProbeSets() {
		CHPFileHeader h = new CHPFileHeader();
		h.setNumProbeSets(10);
		assertEquals(h.getNumProbeSets(), 10);
	}

	/**
	 * Test of getAssayType method, of class affymetrix.gcos.chp.CHPFileHeader.
	 */
	public void testAssayType() {
		CHPFileHeader h = new CHPFileHeader();
		h.setAssayType(CHPFileHeader.EXPRESSION_ASSAY);
		assertEquals(h.getAssayType(), CHPFileHeader.EXPRESSION_ASSAY);
	}

	/**
	 * Test of getChipType method, of class affymetrix.gcos.chp.CHPFileHeader.
	 */
	public void testChipType() {
		CHPFileHeader h = new CHPFileHeader();
		h.setChipType("Test3");
		assertEquals(h.getChipType(), "Test3");
	}

	/**
	 * Test of getAlgName method, of class affymetrix.gcos.chp.CHPFileHeader.
	 */
	public void testAlgName() {
		CHPFileHeader h = new CHPFileHeader();
		h.setAlgName("A1");
		assertEquals(h.getAlgName(), "A1");
	}

	/**
	 * Test of getAlgVersion method, of class affymetrix.gcos.chp.CHPFileHeader.
	 */
	public void testAlgVersion() {
		CHPFileHeader h = new CHPFileHeader();
		h.setAlgVersion("v12");
		assertEquals(h.getAlgVersion(), "v12");
	}

	/**
	 * Test of getParentCellFile method, of class affymetrix.gcos.chp.CHPFileHeader.
	 */
	public void testParentCellFile() {
		CHPFileHeader h = new CHPFileHeader();
		h.setParentCellFile("C:\\test.cel");
		assertEquals(h.getParentCellFile(), "C:\\test.cel");

	}

	/**
	 * Test of getProgID method, of class affymetrix.gcos.chp.CHPFileHeader.
	 */
	public void testProgID() {
		CHPFileHeader h = new CHPFileHeader();
		h.setProgID("id");
		assertEquals(h.getProgID(), "id");

	}

	/**
	 * Test of getAlgorithmParameters method, of class affymetrix.gcos.chp.CHPFileHeader.
	 */
	public void testAlgorithmParameters() {
		CHPFileHeader h = new CHPFileHeader();
		List<TagValuePair> params = new ArrayList<TagValuePair>();
		TagValuePair p = new TagValuePair();
		p.setTag("t1");
		p.setValue("v1");
		params.add(p);
		p = new TagValuePair();
		p.setTag("t2");
		p.setValue("v2");
		params.add(p);
		h.setAlgorithmParameters(params);

		List<TagValuePair> params2 = h.getAlgorithmParameters();
		p = params2.get(0);
		assertEquals(p.getTag(), "t1");
		assertEquals(p.getValue(), "v1");
		p = params2.get(1);
		assertEquals(p.getTag(), "t2");
		assertEquals(p.getValue(), "v2");
	}

	/**
	 * Test of getSummaryParameters method, of class affymetrix.gcos.chp.CHPFileHeader.
	 */
	public void testSummaryParameters() {
		CHPFileHeader h = new CHPFileHeader();
		List<TagValuePair> params = new ArrayList<TagValuePair>();
		TagValuePair p = new TagValuePair();
		p.setTag("t1");
		p.setValue("v1");
		params.add(p);
		p = new TagValuePair();
		p.setTag("t2");
		p.setValue("v2");
		params.add(p);
		h.setSummaryParameters(params);

		List<TagValuePair> params2 = h.getSummaryParameters();
		p = params2.get(0);
		assertEquals(p.getTag(), "t1");
		assertEquals(p.getValue(), "v1");
		p = params2.get(1);
		assertEquals(p.getTag(), "t2");
		assertEquals(p.getValue(), "v2");

	}

	/**
	 * Test of getBackgroundZoneInfo method, of class affymetrix.gcos.chp.CHPFileHeader.
	 */
	public void testBackgroundZoneInfo() {
		CHPFileHeader h = new CHPFileHeader();
		BackgroundZoneInfo z = new BackgroundZoneInfo();
		z.setSmoothFactor(1.0f);
		h.setBackgroundZoneInfo(z);
		BackgroundZoneInfo z2 = h.getBackgroundZoneInfo();
		assertEquals(z2.getSmoothFactor(), 1.0f, 0.000001f);

	}

	/**
	 * Test of get summary parameters.
	 */
	public void testGetSummaryParameter() {
		CHPFileHeader h = new CHPFileHeader();
		TagValuePair p = new TagValuePair();
		List<TagValuePair> params = new ArrayList<TagValuePair>();
		p.setTag("tag1");
		p.setValue("value1");
		params.add(p);
		h.setSummaryParameters(params);
		assertEquals(h.getSummaryParameter("tag1"), "value1");
		assertEquals(h.getSummaryParameter("tag2"), null);
	}

	/**
	 * Test of get alg parameters.
	 */
	public void testGetAlgorithmParameters() {
		CHPFileHeader h = new CHPFileHeader();
		TagValuePair p = new TagValuePair();
		List<TagValuePair> params = new ArrayList<TagValuePair>();
		p.setTag("tag1");
		p.setValue("value1");
		params.add(p);
		h.setAlgorithmParameters(params);
		assertEquals(h.getAlgorithmParameter("tag1"), "value1");
		assertEquals(h.getAlgorithmParameter("tag2"), null);
	}

	/**
	 * Test of get bg zone.
	 */
	public void testGetBackgroundZone() {
		// getBackgroundZone

		CHPFileHeader h = new CHPFileHeader();

		BackgroundZoneInfo z = new BackgroundZoneInfo();
		assertEquals(z.getNumberZones(), 0);
		z.addZone(new BackgroundZoneType(1.0f, 2.0f, 3.0f));
		z.addZone(new BackgroundZoneType(11.0f, 12.0f, 13.0f));
		z.addZone(new BackgroundZoneType(21.0f, 22.0f, 23.0f));
		assertEquals(z.getNumberZones(), 3);

		h.setBackgroundZoneInfo(z);

		BackgroundZoneType zone = h.getBackgroundZone(11, 12);

		assertEquals(zone.getCenterX(), 11.0f, 0.00001f);
		assertEquals(zone.getCenterY(), 12.0f, 0.00001f);
		assertEquals(zone.getBackground(), 13.0f, 0.00001f);

	}
}
