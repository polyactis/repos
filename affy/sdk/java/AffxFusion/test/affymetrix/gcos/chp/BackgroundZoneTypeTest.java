/*
 * BackgroundZoneTypeTest.java
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
public class BackgroundZoneTypeTest extends TestCase {

	public BackgroundZoneTypeTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(BackgroundZoneTypeTest.class);

		return suite;
	}

	/**
	 * Test of getCenterX method, of class affymetrix.gcos.chp.BackgroundZoneType.
	 */
	public void testCenterX() {
		BackgroundZoneType z = new BackgroundZoneType();
		z.setCenterX(10.0f);
		assertEquals(z.getCenterX(), 10.0f, 0.000001f);
	}

	/**
	 * Test of getCenterY method, of class affymetrix.gcos.chp.BackgroundZoneType.
	 */
	public void testCenterY() {
		BackgroundZoneType z = new BackgroundZoneType();
		z.setCenterY(10.0f);
		assertEquals(z.getCenterY(), 10.0f, 0.000001f);
	}

	/**
	 * Test of getBackground method, of class affymetrix.gcos.chp.BackgroundZoneType.
	 */
	public void testBackground() {
		BackgroundZoneType z = new BackgroundZoneType();
		z.setBackground(10.0f);
		assertEquals(z.getBackground(), 10.0f, 0.000001f);
	}

	/**
	 * Test of constructor.
	 */
	public void testConstructor() {
		BackgroundZoneType z1 = new BackgroundZoneType();
		z1.setCenterX(5.0f);
		z1.setCenterY(10.0f);
		z1.setBackground(1.0f);
		BackgroundZoneType z2 = new BackgroundZoneType(z1);
		z1.setCenterX(0.0f);
		z1.setCenterY(0.0f);
		z1.setBackground(0.0f);

		assertEquals(z2.getCenterX(), 5.0f, 0.000001f);
		assertEquals(z2.getCenterY(), 10.0f, 0.000001f);
		assertEquals(z2.getBackground(), 1.0f, 0.000001f);
	}

}
