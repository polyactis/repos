/*
 * BackgroundZoneInfoTest.java
 * JUnit based test
 *
 * Created on October 12, 2005, 8:15 PM
 */

package affymetrix.gcos.chp;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * 
 * @author ljevon
 */
public class BackgroundZoneInfoTest extends TestCase {

	public BackgroundZoneInfoTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(BackgroundZoneInfoTest.class);

		return suite;
	}

	/**
	 * Test of getNumberZones method, of class affymetrix.gcos.chp.BackgroundZoneInfo.
	 */
	public void testGetNumberZones() {
		BackgroundZoneInfo z = new BackgroundZoneInfo();
		assertEquals(z.getNumberZones(), 0);
		z.addZone(new BackgroundZoneType(1.0f, 2.0f, 3.0f));
		assertEquals(z.getNumberZones(), 1);
		z.clear();
		assertEquals(z.getNumberZones(), 0);
	}

	/**
	 * Test of addZone method, of class affymetrix.gcos.chp.BackgroundZoneInfo.
	 */
	public void testZone() {
		BackgroundZoneInfo z = new BackgroundZoneInfo();
		assertEquals(z.getNumberZones(), 0);
		z.addZone(new BackgroundZoneType(1.0f, 2.0f, 3.0f));
		z.addZone(new BackgroundZoneType(11.0f, 12.0f, 13.0f));
		z.addZone(new BackgroundZoneType(21.0f, 22.0f, 23.0f));
		assertEquals(z.getNumberZones(), 3);
		BackgroundZoneType zone;
		zone = z.getZone(0);
		assertEquals(zone.getCenterX(), 1.0f, 0.00001f);
		assertEquals(zone.getCenterY(), 2.0f, 0.00001f);
		assertEquals(zone.getBackground(), 3.0f, 0.00001f);
		zone = z.getZone(1);
		assertEquals(zone.getCenterX(), 11.0f, 0.00001f);
		assertEquals(zone.getCenterY(), 12.0f, 0.00001f);
		assertEquals(zone.getBackground(), 13.0f, 0.00001f);
		zone = z.getZone(2);
		assertEquals(zone.getCenterX(), 21.0f, 0.00001f);
		assertEquals(zone.getCenterY(), 22.0f, 0.00001f);
		assertEquals(zone.getBackground(), 23.0f, 0.00001f);
	}

	/**
	 * Test of getSmoothFactor method, of class affymetrix.gcos.chp.BackgroundZoneInfo.
	 */
	public void testSmoothFactor() {
		BackgroundZoneInfo z = new BackgroundZoneInfo();
		z.setSmoothFactor(1.0f);
		assertEquals(z.getSmoothFactor(), 1.0f, 0.000001f);
	}

}
