/*
 * FusionCELFileEntryTypeTest.java
 * JUnit based test
 *
 * Created on October 16, 2005, 8:04 PM
 */

package affymetrix.fusion.cel;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * 
 * @author ljevon
 */
public class FusionCELFileEntryTypeTest extends TestCase {

	public FusionCELFileEntryTypeTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(FusionCELFileEntryTypeTest.class);

		return suite;
	}

	/**
	 * Test of getIntensity method, of class affymetrix.gcos.cel.CELFileEntryType.
	 */
	public void testIntensity() {
		FusionCELFileEntryType c = new FusionCELFileEntryType();
		c.setIntensity(10.0f);
		assertEquals(c.getIntensity(), 10.0f, 0.00000001f);
		assertEquals(c.getStdv(), 0.0f, 0.00000001f);
		assertEquals(c.getPixels(), 0);
	}

	/**
	 * Test of getStdv method, of class affymetrix.gcos.cel.CELFileEntryType.
	 */
	public void testStdv() {
		FusionCELFileEntryType c = new FusionCELFileEntryType();
		c.setStdv(10.0f);
		assertEquals(c.getIntensity(), 0.0f, 0.00000001f);
		assertEquals(c.getStdv(), 10.0f, 0.00000001f);
		assertEquals(c.getPixels(), 0);
	}

	/**
	 * Test of getPixels method, of class affymetrix.gcos.cel.CELFileEntryType.
	 */
	public void testPixels() {
		FusionCELFileEntryType c = new FusionCELFileEntryType();
		c.setPixels((short)10);
		assertEquals(c.getIntensity(), 0.0f, 0.00000001f);
		assertEquals(c.getStdv(), 0.0f, 0.00000001f);
		assertEquals(c.getPixels(), 10);
	}

}
